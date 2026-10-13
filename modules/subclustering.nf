process SUBCLUSTERING {
    tag "Subsetting & Reclustering: ${meta.id}"
    publishDir "${params.outdir}/08_subclustering", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_annotated)
    path interest_csv
    path custom_marker_file
    path marker_file
    val(annotated)
    val(tissue_name)
    
    output:
    tuple val(meta), path("${meta.id}.*_adata_clustered.h5ad"),  emit: adata
    tuple val(meta), path("${meta.id}.*_cluster_sizes.csv"),     emit: cluster_sizes
    tuple val(meta), path("${meta.id}.*_hvg_plot.pdf"),          emit: hvg_plot
    tuple val(meta), path("${meta.id}.*_pca_variance.pdf"),      emit: pca_plot
    tuple val(meta), path("${meta.id}.*_sctype_scores.csv"),     emit: sctype_scores
    tuple val(meta), path("${meta.id}.*.pdf"),          emit: visualization_plots
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import anndata
    import re
    import matplotlib.pyplot as plt
    import seaborn as sns
    import sys
    import os
    # Dynamically add the directory containing sctypepy to the Python Path
    sys.path.append(os.path.join("${workflow.projectDir}", "bin"))
    from sctypepy.sctype_py import gene_sets_prepare, sctype_score, process_cluster

    sc.settings.verbosity = 3

    # Natively read configuration constants as explicit numeric types with safe Nextflow fallbacks
    py_sub_n_pcs      = int("${params.sub_n_pcs ?: 20}")
    py_sub_neighbors  = int("${params.sub_neighbors ?: 10}")
    py_sub_pcs_umap   = int("${params.sub_pcs_umap ?: 15}")
    py_sub_resolution = float("${params.sub_resolution ?: 0.5}")
    py_tissue_name    = "${params.tissue_name ?: 'Eye'}".strip()

    # ==========================================
    # 1. READ ANNDATA FIRST
    # ==========================================
    print("Loading dimensionality reduction data for ${meta.id}...")
    adata = sc.read_h5ad("${adata_annotated}")

    # ==========================================
    # 2. READ CSV & EXTRACT CELLS OF INTEREST
    # ==========================================
    print("Reading target cells of interest from CSV...")
    df_interest = pd.read_csv('${interest_csv}')
    df_interest.columns = df_interest.columns.str.strip()
    
    clean_id = "${meta.id}".replace("_adata", "").strip()
    id_col_name = df_interest.columns[0]
    df_sample = df_interest[df_interest[id_col_name].astype(str).str.strip() == clean_id]
    
    if df_sample.empty:
        raise ValueError(f"Error: Sample ID '{clean_id}' derived from '${meta.id}' not found in CSV column '{id_col_name}'!")

    target_column_in_obs = df_sample.iloc[0, 1].strip()          
    raw_identities = str(df_sample.iloc[0, 2]).strip()        
    allowed_identities = [item.strip() for item in raw_identities.split(',')]
    
    identities_str = "_".join([str(x) for x in allowed_identities])
    identities_safe = re.sub(r'[^a-zA-Z0-9_\\-]', '', identities_str.replace(' ', '_'))
    file_prefix = f"${meta.id}.{identities_safe}"
    
    if target_column_in_obs not in adata.obs.columns:
        raise ValueError(f"Column '{target_column_in_obs}' not found in adata.obs metadata!")

    adata = adata[adata.obs[target_column_in_obs].isin(allowed_identities)].copy()
    if adata.n_obs == 0:
        raise ValueError("Error: Zero cells left after filtering down with interest_csv!")

    # ==========================================
    # 3. RERUN PIPELINE (HVG -> DIM RED -> CLUSTER)
    # ==========================================
    print("Identify HVG...")
    adata.X = adata.layers["scran_normalization"].copy()
    sc.pp.highly_variable_genes(adata, layer="scran_normalization")

    plt.figure(figsize=(6, 5))
    ax = sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=15)
    ax.set_xlim(None, 1.5)
    ax.set_ylim(None, 3)
    plt.savefig(f'{file_prefix}_hvg_plot.pdf', format="pdf", bbox_inches="tight", dpi=300)
    plt.close()

    # PCA
    print("Computing PCA...")
    adata.var["highly_variable"] = adata.var["highly_deviant"].copy()
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable", n_comps=py_sub_n_pcs)

    # Plot variance explained
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=py_sub_n_pcs, show=False)
    fig = plt.gcf()
    fig.suptitle(f"PCA Variance Ratio - {clean_id}", y=1.02, fontsize=12)
    plt.savefig(f'{file_prefix}_pca_variance.pdf', format="pdf", bbox_inches="tight", dpi=300)
    plt.close()

    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=py_sub_neighbors, n_pcs=py_sub_pcs_umap)

    print("Computing UMAP and t-SNE...")
    sc.tl.umap(adata)
    sc.tl.tsne(adata, use_rep='X_pca')

    print("Performing Leiden clustering...")
    sc.tl.leiden(adata, flavor="igraph", resolution=py_sub_resolution, n_iterations=2)

    if 'dendrogram_leiden' in adata.uns:
        print("Clearing inherited parent dendrogram mapping entries...")
        del adata.uns['dendrogram_leiden']

    # Compute a clean, matching sub-dendrogram mapping sequence for sub-clusters
    print("Computing a fresh subcluster hierarchical dendrogram linkage map...")
    sc.tl.dendrogram(adata, groupby='leiden', use_rep='X_pca')

    cluster_df = adata.obs['leiden'].value_counts().reset_index()
    cluster_df.columns = ['cluster', 'n_cells']
    cluster_df = cluster_df.sort_values('cluster')
    cluster_df.to_csv(f'{file_prefix}_cluster_sizes.csv', index=False)

    # ==========================================
    # 4. CELL TYPE ANNOTATION
    # ==========================================
    print("Annotating cell types...")
    gs_list = gene_sets_prepare('${marker_file}', cell_type=py_tissue_name)
    # Scale the data before scoring
    sc.pp.scale(adata, max_value=10)
    scaled_data = pd.DataFrame(adata.X)
    # change column indexes
    scaled_data.columns = adata.var_names
    # Change the row indexes
    scaled_data.index = adata.obs_names
    scRNAseqData = scaled_data.T
    es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])
    unique_clusters = adata.obs['leiden'].unique()

    cL_results = pd.concat([process_cluster(cluster, adata, es_max, 'leiden') for cluster in unique_clusters])
    sctype_scores = cL_results.sort_values('scores', ascending=False).drop_duplicates(subset=['cluster']).sort_values('type').reset_index(drop=True)
    sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

    adata.obs['sctype_classification'] = ""
    for cluster in sctype_scores['cluster'].unique():
        cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
        cl_type_value = cl_type['type'].iloc[0]
        adata.obs.loc[adata.obs['leiden'] == cluster, 'sctype_classification'] = cl_type_value

    output_h5ad = f'{file_prefix}_adata_clustered.h5ad'
    adata.write(output_h5ad)
    sctype_scores.to_csv(f'{file_prefix}_sctype_scores.csv', index=False)
    print(f"Subclustering and annotation complete for ${meta.id}")


    from run_visualization import run_visualization
    
    raw_custom_file = "${custom_marker_file}".strip()
    primary_marker = "${marker_file}".strip()
    
    # 🟢 FIXED: Never fall back to a binary .xlsx file for text parsing loops
    # If the user-supplied custom CSV exists and isn't a Nextflow empty list mock, use it.
    if raw_custom_file and raw_custom_file != "[:]" and os.path.exists(raw_custom_file):
        active_marker_file = raw_custom_file
    elif primary_marker.endswith('.csv') and os.path.exists(primary_marker):
        active_marker_file = primary_marker
    else:
        # Fallback to an empty dictionary state or a known valid CSV layout to prevent binary crash hooks
        print("Warning: No valid CSV marker file found for visualization. Generating placeholder.")
        active_marker_file = "empty_markers.csv"
        with open(active_marker_file, 'w') as f:
            f.write("Placeholder,NoMarkersFound\\n")

    run_visualization(
        adata_path         = output_h5ad,
        sample_id          = file_prefix,
        marker_genes_file  = active_marker_file,
        annotated_str      = "${annotated}"
    )

    """
}
