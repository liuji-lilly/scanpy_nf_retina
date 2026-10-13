process MARKER_GENES {
    tag "Marker gene identification: ${meta.id}"
    publishDir "${params.outdir}/05_markers", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_clustered)
    
    output:
    tuple val(meta), path("${meta.id}_adata_with_markers.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_marker_genes.csv"), emit: markers
    tuple val(meta), path("*.png"), emit: plots
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white')
    
    # Load data using the dynamic metadata dictionary pointer
    print(f"Loading clustered data for ${meta.id}...")
    adata = sc.read_h5ad('${adata_clustered}')
    
    # Find marker genes
    print("Identifying marker genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Extract top markers
    n_genes = 20
    markers_list = []
    
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster, key='rank_genes_groups')
        cluster_markers['cluster'] = cluster
        cluster_markers = cluster_markers.head(n_genes)
        markers_list.append(cluster_markers)
    
    markers_df = pd.concat(markers_list)
    markers_df.to_csv('${meta.id}_marker_genes.csv', index=False)
    
    # ENLARGE TEXT GLOBALLY: Changes x-tick labels, y-tick labels, and titles
    matplotlib.rcParams['xtick.labelsize'] = 12
    matplotlib.rcParams['ytick.labelsize'] = 12
    matplotlib.rcParams['axes.labelsize'] = 14

    # =========================================================================
    # DYNAMIC SIZING FOR AUTOMATED PLOTS (Top 10 genes per cluster)
    # =========================================================================
    num_clusters = len(adata.obs['leiden'].cat.categories)
    approx_top_genes = num_clusters * 10
    auto_width = max(10, approx_top_genes * 0.35 + 4) 
    auto_height = max(6, num_clusters * 0.4 + 2)

    # Plot heatmap of top markers
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='leiden', 
                                     show_gene_labels=True, show=False, 
                                     figsize=(auto_width, auto_height))
    plt.savefig('${meta.id}_marker_genes_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Plot dotplot of top markers
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, groupby='leiden', 
                                     figsize=(auto_width, auto_height))
    plt.savefig('${meta.id}_marker_genes_dotplot.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save with unique dynamic naming
    adata.write('${meta.id}_adata_with_markers.h5ad')
    print(f"Marker gene identification complete for ${meta.id}")
    """
}
