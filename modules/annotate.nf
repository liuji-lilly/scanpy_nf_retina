process ANNOTATE_CELL_TYPES {
    tag "Cell type annotation: ${meta.id}"
    publishDir "${params.outdir}/06_annotation", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_with_markers)
    path(marker_file)       
    val(tissue_name)        
    val(species)
    output:
    tuple val(meta), path("${meta.id}_adata_annotated.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_sctype_scores.csv"), emit: sctype_scores
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    import os

    # Dynamically add the directory containing sctypepy to the Python Path
    sys.path.append(os.path.join("${workflow.projectDir}", "bin"))

    import scanpy as sc
    import pandas as pd
    from sctypepy.sctype_py import gene_sets_prepare, sctype_score, process_cluster
    sc.settings.verbosity = 3
    
    # Load data
    print(f"Loading data for sample: ${meta.id}...")
    adata = sc.read_10x_mtx('${adata_with_markers}') if os.path.isdir('${adata_with_markers}') else sc.read_h5ad('${adata_with_markers}')
    
    # Track the requested species string
    target_species = '${species}'.strip().lower()
    print(f"Processing dataset using casing rules for species: {target_species}")

    # 1. Normalize AnnData dataset index symbols based on target species
    if target_species == 'human':
        adata.var_names = [gene.upper() for gene in adata.var_names]
    elif target_species == 'mouse':
        adata.var_names = [gene.capitalize() for gene in adata.var_names]
    else:
        print(f"WARNING: Unknown species '{target_species}'. Retaining original matrix casing.")
    
    # Prepare gene sets from custom markers file
    gs_list = gene_sets_prepare('${marker_file}', cell_type='${tissue_name}')

    # 2. Normalize ScType parsed markers list to guarantee structural alignment matches
    for list_key in ['gs_positive', 'gs_negative']:
        normalized_lists = []
        for marker_sublist in gs_list[list_key]:
            if target_species == 'human':
                normalized_lists.append([gene.upper() for gene in marker_sublist])
            elif target_species == 'mouse':
                normalized_lists.append([gene.capitalize() for gene in marker_sublist])
            else:
                normalized_lists.append(marker_sublist)
        gs_list[list_key] = normalized_lists

    # Scale the data before scoring
    sc.pp.scale(adata, max_value=10)
    scaled_data = pd.DataFrame(adata.X)
    # change column indexes
    scaled_data.columns = adata.var_names
    # Change the row indexes
    scaled_data.index = adata.obs_names
    scRNAseqData = scaled_data.T

    # Calculate ScType scores for every individual cell
    es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])
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

    # Save with unique dynamic naming
    adata.write('${meta.id}_adata_annotated.h5ad')
    sctype_scores.to_csv('${meta.id}_sctype_scores.csv', index=False)
    print(f"Cell type annotation complete for ${meta.id}")
    """
}