process VISUALIZATION {
    tag "Generate visualizations: ${meta.id}"
    publishDir "${params.outdir}/07_visualization", mode: 'copy'
    
    input:
    tuple val(meta), path(adata)
    val(annotated)
    path marker_genes_file
    
    output:
    tuple val(meta), path("*.pdf"), emit: plots
    
    script:
    """
    #!/usr/bin/env python3
    import csv
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(facecolor='white')

    print("Loading data for visualization on sample: ${meta.id}...")
    adata = sc.read_h5ad('${adata}')
    num_clusters = len(adata.obs['leiden'].unique())

    # =========================================================================
    # GLOBAL PALETTE CONFIGURATION
    # =========================================================================
    lillycolor = ["#0F3A85","#8A969E","#99BFE5","#E1251B", "#212121", "#FFC709","#511207", "#F58E7D", "#E4EBF1", "#FDE8E5", "#144B2D", "#C6DCD8"]
    pal1 = ["#8B633DFF", "#8A999DFF", "#C4B57CFF", "#CB4E33FF", "#F0D1A9FF", "#1B3F8AFF", "#68A1D3FF", "#AECCDFFF", "#D18860FF", "#30588BFF", "#C2944FFF", "#F9D073FF"]
    pal2 = ["#C70E7BFF", "#FC6882FF", "#007BC3FF", "#54BCD1FF", "#EF7C12FF", "#F4B95AFF", "#009F3FFF", "#8FDA04FF", "#AF6125FF", "#F4E3C7FF", "#B25D91FF", "#EFC7E6FF", "#EF7C12FF", "#F4B95AFF"]
    pal3 = ["#402733FF", "#8F6052FF", "#2F6349FF", "#AE2B47FF", "#7BC17DFF", "#68D8F2FF", "#7982F7FF", "#E03E4AFF", "#FFF759FF", "#ECA987FF", "#705549FF", "#9F4739FF"]
    pal4 = ["#3B2319FF", "#80521CFF", "#D29C44FF", "#EBC174FF", "#EDE2CCFF", "#7EC5F4FF", "#4585B7FF", "#225E92FF", "#183571FF", "#43429BFF", "#5E65BEFF"]
    pal5 = ["#FF0066FF", "#328C97FF", "#D1AAC2FF", "#A5506DFF", "#B3E0BFFF", "#2A9D3DFF", "#EDF181FF", "#DB7003FF", "#FBA600FF", "#F8C1A6FF", "#A30000FF", "#FF3200FF", "#011A51FF", "#97D1D9FF", "#916C37FF"]
    
    customcolor = lillycolor + pal1 + pal2 + pal3 + pal4 + pal5
    custom_cmap = ListedColormap(customcolor)

     # Clean and order categories so Scanpy assigns colors in your precise order
    if 'leiden' in adata.obs.columns:
        # Enforce natural numerical sorting (0, 1, 2, ... 10, 11) for cluster mapping
        adata.obs['leiden'] = pd.Categorical(adata.obs['leiden'], categories=sorted(adata.obs['leiden'].unique(), key=int), ordered=True)
        adata.uns['leiden_colors'] = customcolor[:len(adata.obs['leiden'].cat.categories)]

    target_obs = None
    if 'sctype_classification' in adata.obs.columns:
        target_obs = 'sctype_classification'
    elif 'cell_type' in adata.obs.columns:
        target_obs = 'cell_type'

    if target_obs:
        # Enforce consistent sorting for cell type labels
        adata.obs[target_obs] = pd.Categorical(adata.obs[target_obs], categories=sorted(adata.obs[target_obs].unique()), ordered=True)
        adata.uns[f'{target_obs}_colors'] = customcolor[:len(adata.obs[target_obs].cat.categories)]

    # =========================================================================
    # 1. CORE EMBEDDING PLOTS
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), layout='constrained')
    sc.pl.umap(adata, color=target_obs, title="UMAP - Cell Types (${meta.id})", show=False, legend_loc='right margin',ax=axes[0])
    sc.pl.umap(adata, color=target_obs, title="UMAP - Cell Types (${meta.id})", show=False, legend_loc='on data',legend_fontweight='bold',ax=axes[1])
    plt.savefig('${meta.id}_umap_cell_types.pdf', bbox_inches='tight')
    plt.close()

    if "${annotated}" == "true" and target_obs:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5), layout='constrained')
        sc.pl.umap(adata, color=target_obs, title="UMAP - Cell Types (${meta.id})", show=False, legend_loc='right margin',ax=axes[0])
        sc.pl.umap(adata, color=target_obs, title="UMAP - Cell Types (${meta.id})", show=False, legend_loc='on data',legend_fontweight='bold',ax=axes[1])
        plt.savefig('${meta.id}_umap_cell_types.pdf', bbox_inches='tight')
        plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), layout='constrained')
    sc.pl.pca(adata, color='leiden', title="PCA - Clusters (${meta.id})", show=False, legend_loc='right margin',annotate_var_explained=True,ax=axes[0])
    sc.pl.pca(adata, color='leiden', title=f"PCA - Clusters (${meta.id})", show=False, legend_loc='on data',legend_fontweight='bold', annotate_var_explained=True,ax=axes[1])
    plt.savefig('${meta.id}_pca_clusters.pdf', bbox_inches='tight')
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), layout='constrained')
    sc.pl.tsne(adata, color='leiden', title="t-SNE - Clusters (${meta.id})", show=False, legend_loc='right margin',ax=axes[0])
    sc.pl.tsne(adata, color='leiden', title="t-SNE - Clusters (${meta.id})", show=False, legend_loc='on data',legend_fontweight='bold',ax=axes[1])
    plt.savefig('${meta.id}_tsne_clusters.pdf', bbox_inches='tight')
    plt.close()

    # =========================================================================
    # 2. QC METRICS CONSOLIDATED GRIDS
    # =========================================================================
    qc_cols = ["leiden", "total_counts", "pct_counts_mt", "n_genes_by_counts"]
    valid_qc = [c for c in qc_cols if c in adata.obs.columns]
    if "predicted_doublet" in adata.obs.columns: valid_qc.append("predicted_doublet")
    if "doublet_score" in adata.obs.columns: valid_qc.append("doublet_score")

    # Continuous parameters default to sequential maps, categorical values obey customcolor automatically
    sc.pl.umap(adata, color=valid_qc, wspace=0.4, ncols=2, show=False)
    plt.savefig('${meta.id}_umap_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    sc.pl.pca(adata, color=valid_qc, wspace=0.4, ncols=2, show=False)
    plt.savefig('${meta.id}_pca_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    sc.pl.tsne(adata, color=valid_qc, wspace=0.4, ncols=2, show=False)
    plt.savefig('${meta.id}_tsne_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    # =========================================================================
    # 3. RETINA MARKER GENE DIGESTS
    # =========================================================================
    print("Loaded marker dictionary:")
    marker_genes_dict = {}

    with open('${marker_genes_file}', mode="r", encoding="utf-8") as file:
        reader = csv.reader(file)
        for row in reader:
            if row:
                cell_type = row[0]
                genes = row[1:]
                marker_genes_dict[cell_type] = genes

    available_markers = {}
    flat_gene_list = []
    for cell_type, genes in marker_genes_dict.items():
        available = [g for g in genes if g in adata.var_names]
        if available:
            available_markers[cell_type] = available
            for g in available:
                if g not in flat_gene_list:
                    flat_gene_list.append(g)

    if available_markers:
        custom_width = max(6, len(flat_gene_list) * 0.45 + 3)
        custom_height = max(4, num_clusters * 0.4 + 2)

        # Dotplots and Stacked Violins capture the registered uns color arrays gracefully
        sc.pl.dotplot(adata, available_markers, "leiden", dendrogram=True, show=False, figsize=(custom_width, custom_height))
        plt.savefig('${meta.id}_marker_selectedmarkers_dotplot.pdf', bbox_inches='tight')
        plt.close()

        sc.pl.stacked_violin(adata, available_markers, groupby="leiden", swap_axes=False, dendrogram=True, show=False, figsize=(custom_width, custom_height))
        plt.savefig('${meta.id}_marker_selectedmarkers_stacked_violin.pdf', bbox_inches='tight')
        plt.close()

        sc.pl.umap(adata, color=flat_gene_list, s=5, frameon=False, ncols=4, vmin=0, vmax="p99", show=False)
        plt.savefig('${meta.id}_marker_custom_umap_expression.pdf', bbox_inches='tight')
        plt.close()

        combined_features = ['leiden'] + flat_gene_list
        sc.pl.pca(adata, color=combined_features, s=6, frameon=False, ncols=4, wspace=0.3, vmin=0, vmax="p99", show=False)
        plt.savefig('${meta.id}_combined_pca_clusters_and_markers.pdf', bbox_inches='tight')
        plt.close()

    # =========================================================================
    # 4. CELL TYPE / CLUSTER DISTRIBUTION (STACKED BAR PLOTS)
    # =========================================================================
    if target_obs:
        print("Generating stacked bar plots for cell type distributions...")
        
        # 1. Global Sample Composition (Single Stacked Bar)
        counts = adata.obs[target_obs].value_counts()
        proportions = (counts / counts.sum() * 100).to_frame().T
        
        # Get the exactly aligned discrete list of custom colors for the present cell types
        num_cell_types = proportions.shape[1]
        global_colors = customcolor[:num_cell_types]
        
        fig, ax = plt.subplots(figsize=(5, 8))
        # FIX: Changed cmap=custom_cmap to color=global_colors to preserve discrete Lilly colors
        proportions.plot(kind='bar', stacked=True, ax=ax, edgecolor='none', color=global_colors)
        
        ax.grid(False)
        ax.set_ylabel('Percentage of Sample (%)')
        ax.set_xlabel('${meta.id}')
        ax.set_xticklabels([]) 
        
        plt.legend(title='Cell Types', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig('${meta.id}_global_celltype_composition.pdf', bbox_inches='tight')
        plt.close()

        # 2. Representation of Leiden Clusters within each Cell Type Category
        cluster_in_ct = pd.crosstab(adata.obs[target_obs], adata.obs['leiden'], normalize='index') * 100
        
        # Get the exactly aligned discrete list of custom colors for the present Leiden clusters
        num_clusters_present = cluster_in_ct.shape[1]
        cluster_colors = customcolor[:num_clusters_present]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        # FIX: Changed cmap=custom_cmap to color=cluster_colors to preserve discrete Lilly colors
        cluster_in_ct.plot(kind='bar', stacked=True, ax=ax, edgecolor='none', color=cluster_colors)
        
        ax.grid(False)
        ax.set_title("Cluster Composition per Cell Type (${meta.id})")
        ax.set_ylabel('Percentage (%)')
        ax.set_xlabel('Cell Type Classification')
        plt.xticks(rotation=30, ha='right')
        
        plt.legend(title='Leiden Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig('${meta.id}_cluster_composition_per_celltype.pdf', bbox_inches='tight')
        plt.close()    
        
    else:
        print("Skipping distribution plots: target metadata column missing from data matrix obs.")
    """
}
