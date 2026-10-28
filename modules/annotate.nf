process ANNOTATE_CELL_TYPES {
    tag "Cell type annotation: ${meta.id}"
    publishDir "${params.outdir}/06_annotation", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_with_markers)
    path(marker_file)       
    val(tissue_name)        
    output:
    tuple val(meta), path("${meta.id}_adata_annotated.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_sctype_scores.csv"), emit: sctype_scores
    tuple val(meta), path("${meta.id}_cell_type_summary.csv"), emit: summary_csv

    script:
    """
    #!/usr/bin/env python3
    import sys
    import os

    # Dynamically add the directory containing sctypepy to the Python Path
    sys.path.append(os.path.join("${workflow.projectDir}", "bin"))

    import scanpy as sc
    import pandas as pd
    import numpy as np
    from sctypepy.sctype_py import gene_sets_prepare, sctype_score, process_cluster
    sc.settings.verbosity = 3
    import scipy.sparse as sp
    from sklearn.preprocessing import scale
    from collections import defaultdict
    
    def gene_sets_prepare_native(path_to_db_file, cell_type):
        cell_markers = pd.read_excel(path_to_db_file)
        cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
        
        for col in ['geneSymbolmore1', 'geneSymbolmore2']:
            if col in cell_markers.columns:
                cell_markers[col] = cell_markers[col].fillna('').astype(str).str.replace(" ", "").str.upper()
                
        def clean_gene_list(series):
            valid_rows = [str(val) for val in series if val and str(val).strip() not in ['', 'NAN', 'NONE']]
            if not valid_rows: 
                return []
            all_genes = ','.join(valid_rows).split(',')
            return list(set([g.strip().upper() for g in all_genes if g.strip() != '']))
            
        gs_positive = cell_markers.groupby('cellName')['geneSymbolmore1'].apply(clean_gene_list).to_dict()
        gs_negative = cell_markers.groupby('cellName')['geneSymbolmore2'].apply(clean_gene_list).to_dict()
        return {'gs_positive': gs_positive, 'gs_negative': gs_negative}


    def sctype_score_fixed(scRNAseqData, scaled=True, gs=None, gs2=None, gene_names_to_uppercase=True, *args, **kwargs):
        data_df = scRNAseqData.copy()
        
        if gene_names_to_uppercase:
            data_df.index = data_df.index.str.upper()
            gs = {k: [g.upper() for g in v] for k, v in gs.items()}
            gs2 = {k: [g.upper() for g in v] for k, v in gs2.items()} if gs2 else {}

        marker_stat = defaultdict(int, {gene: sum(gene in genes for genes in gs.values()) for gene in set(gene for genes in gs.values() for gene in genes)})
        marker_sensitivity = pd.DataFrame({'gene_': list(marker_stat.keys()), 'score_marker_sensitivity': list(marker_stat.values())})
        
        min_value = 1
        max_value = len(gs) if len(gs) > 1 else 2
        marker_sensitivity['score_marker_sensitivity'] = 1 - (marker_sensitivity['score_marker_sensitivity'] - min_value) / (max_value - min_value)
        
        names_gs_cp = list(gs.keys())
        names_gs_2_cp = list(gs2.keys()) if gs2 else []
        
        gs_ = {key: [gene for gene in data_df.index if gene in gs[key]] for key in gs}
        gs2_ = {key: [gene for gene in data_df.index if gene in gs2[key]] for key in gs2} if gs2 else {key: [] for key in names_gs_cp}
        
        gs__ = dict(zip(names_gs_cp, gs_.values()))
        gs2__ = dict(zip(names_gs_2_cp, gs2_.values())) if gs2 else {key: [] for key in names_gs_cp}
        
        cell_markers_genes_score = marker_sensitivity[marker_sensitivity['gene_'].isin(set.union(*map(set, gs__.values())))]
        
        # FIX: Recast back into a true Pandas DataFrame if scikit-learn scale turns it into a numpy array
        if not scaled:
            Z_array = scale(data_df.T).T
            Z = pd.DataFrame(Z_array, index=data_df.index, columns=data_df.columns)
        else:
            Z = data_df.copy()
            
        for _, row in cell_markers_genes_score.iterrows():
            if row['gene_'] in Z.index:
                Z.loc[row['gene_']] *= row['score_marker_sensitivity']
                
        marker_genes = list(set().union(*gs__.values(), *gs2__.values()))
        if not marker_genes:
            raise ValueError("❌ Error: Zero marker genes matched your expression dataset matrix rows!")
            
        Z = Z.loc[marker_genes]
        
        # FIX: Corrected Z.shape references to standard indexing
        es = pd.DataFrame(index=gs__.keys(), columns=Z.columns, data=np.zeros((len(gs__), Z.shape[1])))
        
        for gss_, genes in gs__.items():
            if not genes:
                continue
            for j in range(Z.shape[1]):
                gs_z = Z.loc[genes, Z.columns[j]]
                
                has_neg = gs2 and gss_ in gs2__ and len(gs2__[gss_]) > 0
                gz_2 = Z.loc[gs2__[gss_], Z.columns[j]] * -1 if has_neg else pd.Series(dtype=np.float64)
                
                sum_t1 = np.sum(gs_z) / np.sqrt(len(gs_z)) if len(gs_z) > 0 else 0
                sum_t2 = np.sum(gz_2) / np.sqrt(len(gz_2)) if not gz_2.empty and len(gz_2) > 0 else 0
                
                if pd.isna(sum_t2):
                    sum_t2 = 0
                es.loc[gss_, Z.columns[j]] = sum_t1 + sum_t2
                
        return es.dropna(how='all')


    # Load data
    print(f"Loading data for sample: ${meta.id}...")
    adata = sc.read_10x_mtx('${adata_with_markers}') if os.path.isdir('${adata_with_markers}') else sc.read_h5ad('${adata_with_markers}')
    
    # Prepare gene sets from custom markers file
    gs_list = gene_sets_prepare_native('${marker_file}', cell_type='${tissue_name}')

    # Scale the data using scanpy
    sc.pp.scale(adata, max_value=10)

    # Build a clean DataFrame: Transposed (.T) because ScType loops require Genes as ROWS
    dense_X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    scRNAseqData = pd.DataFrame(dense_X, columns=adata.var_names, index=adata.obs_names).T

    # Calculate ScType scores for every individual cell
    es_max = sctype_score_fixed(
        scRNAseqData = scRNAseqData, 
        scaled = True, 
        gs = gs_list['gs_positive'], 
        gs2 = gs_list['gs_negative'],
        gene_names_to_uppercase = True
    )
    unique_clusters = adata.obs['leiden'].unique()

    # Aggregate scores at the cluster level
    cL_results = pd.concat([process_cluster(cluster, adata, es_max, 'leiden') for cluster in unique_clusters])

    # Group by cluster and select the top row based on scores
    sctype_scores = cL_results.sort_values('scores', ascending=False).drop_duplicates(subset=['cluster']).sort_values('type').reset_index(drop=True)

    # Set low-confidence clusters to "Unknown"
    sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

    # Iterate over unique clusters
    adata.obs['sctype_classification'] = ""
    for cluster in sctype_scores['cluster'].unique():
        cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
        cl_type_value = cl_type['type'].iloc[0]
        adata.obs.loc[adata.obs['leiden'] == cluster, 'sctype_classification'] = cl_type_value
    
    # --- Summarize total cell counts and percentage distribution ---
    summary_df = sctype_scores.groupby("type")["ncells"].sum().reset_index()
    total_annotated_cells = summary_df["ncells"].sum()
    summary_df["percentage"] = ((summary_df["ncells"] / total_annotated_cells) * 100).round(2)
    summary_df = summary_df.sort_values(by="ncells", ascending=False).reset_index(drop=True)

    # Save with unique dynamic naming
    adata.write('${meta.id}_adata_annotated.h5ad')
    sctype_scores.to_csv('${meta.id}_sctype_scores.csv', index=False)
    summary_df.to_csv('${meta.id}_cell_type_summary.csv', index=False)
    
    print(f"Cell type annotation complete for ${meta.id}")
    """
}