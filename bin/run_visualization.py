#!/usr/bin/env python3
import csv
import sys
import os
import scanpy as sc
import pandas as pd

import matplotlib
matplotlib.use('Agg')

def run_visualization(adata_path, sample_id, marker_genes_file, annotated_str="false"):
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import scipy.sparse
    import math
    import seaborn as sns
    import numpy as np

    """
    Generates standard diagnostic plots, QC layouts, and marker profiles for scRNA-seq samples.
    """
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(facecolor='white')

    print(f"Loading data for visualization on sample: {sample_id}...")
    adata = sc.read_h5ad(adata_path)
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

    # Clean and order categories so Scanpy assigns colors in your precise order
    if 'leiden' in adata.obs.columns:
        adata.obs['leiden'] = pd.Categorical(adata.obs['leiden'], categories=sorted(adata.obs['leiden'].unique(), key=int), ordered=True)
        adata.uns['leiden_colors'] = customcolor[:len(adata.obs['leiden'].cat.categories)]

    target_obs = None
    if 'sctype_classification' in adata.obs.columns:
        target_obs = 'sctype_classification'
    elif 'cell_type' in adata.obs.columns:
        target_obs = 'cell_type'

    if target_obs:
        adata.obs[target_obs] = pd.Categorical(adata.obs[target_obs], categories=sorted(adata.obs[target_obs].unique()), ordered=True)
        adata.uns[f'{target_obs}_colors'] = customcolor[:len(adata.obs[target_obs].cat.categories)]

    # =========================================================================
    # 1. CORE EMBEDDING PLOTS (With corrected legend widths to prevent overlap)
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(15, 5), layout='constrained')
    plt.subplots_adjust(wspace=0.45)
    sc.pl.umap(adata, color=target_obs, title=f"UMAP - Cell Types ({sample_id})", show=False, legend_loc='right margin', ax=axes[0])
    sc.pl.umap(adata, color=target_obs, title=f"UMAP - Cell Types ({sample_id})", show=False, legend_loc='on data', legend_fontweight='bold', ax=axes[1])
    plt.savefig(f'{sample_id}_umap_cell_types.pdf', bbox_inches='tight')
    plt.close()

    if annotated_str.lower() == "true" and target_obs:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5), layout='constrained')
        plt.subplots_adjust(wspace=0.45)
        sc.pl.umap(adata, color=target_obs, title=f"UMAP - Cell Types ({sample_id})", show=False, legend_loc='right margin', ax=axes[0])
        sc.pl.umap(adata, color=target_obs, title=f"UMAP - Cell Types ({sample_id})", show=False, legend_loc='on data', legend_fontweight='bold', ax=axes[1])
        plt.savefig(f'{sample_id}_umap_cell_types_annotated.pdf', bbox_inches='tight')
        plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(15, 5), layout='constrained')
    plt.subplots_adjust(wspace=0.45)
    sc.pl.pca(adata, color='leiden', title=f"PCA - Clusters ({sample_id})", show=False, legend_loc='right margin', annotate_var_explained=True, ax=axes[0])
    sc.pl.pca(adata, color='leiden', title=f"PCA - Clusters ({sample_id})", show=False, legend_loc='on data', legend_fontweight='bold', annotate_var_explained=False, ax=axes[1])
    plt.savefig(f'{sample_id}_pca_clusters.pdf', bbox_inches='tight')
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(15, 5), layout='constrained')
    plt.subplots_adjust(wspace=0.45)
    sc.pl.tsne(adata, color='leiden', title=f"t-SNE - Clusters ({sample_id})", show=False, legend_loc='right margin', ax=axes[0])
    sc.pl.tsne(adata, color='leiden', title=f"t-SNE - Clusters ({sample_id})", show=False, legend_loc='on data', legend_fontweight='bold', ax=axes[1])
    plt.savefig(f'{sample_id}_tsne_clusters.pdf', bbox_inches='tight')
    plt.close()

    # =========================================================================
    # 2. QC METRICS CONSOLIDATED GRIDS
    # =========================================================================
    qc_cols = ["leiden", "total_counts", "pct_counts_mt", "n_genes_by_counts"]
    valid_qc = [c for c in qc_cols if c in adata.obs.columns]
    if "predicted_doublet" in adata.obs.columns: valid_qc.append("predicted_doublet")
    if "doublet_score" in adata.obs.columns: valid_qc.append("doublet_score")

    sc.pl.umap(adata, color=valid_qc, wspace=0.4, ncols=2, show=False)
    plt.savefig(f'{sample_id}_umap_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    sc.pl.pca(adata, color=valid_qc, wspace=0.4, ncols=2, show=False)
    plt.savefig(f'{sample_id}_pca_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    sc.pl.tsne(adata, color=valid_qc, wspace=0.4, ncols=2, show=False)
    plt.savefig(f'{sample_id}_tsne_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    #Make the companion Violin Plot for valid_qc
    # 1. Force remove any duplicate columns inside adata.obs metadata
    if adata.obs.columns.duplicated().any():
        adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated()]

    # 2. Re-cast to string then category to ensure clean pandas grouping
    adata.obs['leiden'] = adata.obs['leiden'].astype(str).astype('category')

    # 3. Clean valid_qc list to strictly filter out any accidental "leiden" entries
    clean_qc_metrics = [metric for metric in valid_qc if metric != 'leiden']

    # 4. Pull out a clean, isolated dataframe for the requested metrics
    plot_df = sc.get.obs_df(adata, keys=clean_qc_metrics + ['leiden'])
    if plot_df.columns.duplicated().any():
        plot_df = plot_df.loc[:, ~plot_df.columns.duplicated()]

    # 5. Generate the multi-panel layout
    num_qc = len(clean_qc_metrics)
    ncols_qc = 2
    nrows_qc = math.ceil(num_qc / ncols_qc)
    
    # Using layout=None to get rid of the "layout engine incompatible" user warnings
    fig, axes = plt.subplots(
        nrows=nrows_qc, 
        ncols=ncols_qc, 
        figsize=(ncols_qc * 5.5, nrows_qc * 4.0),
        layout=None 
    )
    
    if num_qc == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    # 6. Populate via Seaborn 
    for i, qc_metric in enumerate(clean_qc_metrics):
        sns.violinplot(
            data=plot_df,
            x='leiden',
            y=qc_metric,
            ax=axes[i],
            inner=None,          
            hue='leiden',        
            legend=False
        )
        axes[i].set_title(qc_metric, fontsize=14)
        axes[i].tick_params(axis='x', labelrotation=45)

    # 7. Delete any empty leftover subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # 8. Render and export cleanly
    plt.tight_layout()
    plt.savefig(f'{sample_id}_violin_comprehensive_qc.pdf', bbox_inches='tight')
    plt.close()

    # =========================================================================
    # 3. UNIQUE NON-OVERLAPPING MARKER GENE DIGESTS
    # =========================================================================
    print("Reading marker configurations...")
    marker_genes_dict = {}
    try:
        with open(marker_genes_file, mode="r", encoding="utf-8-sig") as file:
            # Quick check to ensure the file can be read without breaking
            file.readline()
            file.seek(0)
            reader = csv.reader(file)
            rows = list(reader)
    except UnicodeDecodeError:
        print(f"Notice: UTF-8 decoding failed for {marker_genes_file}. Falling back to Windows-1252 (Excel) encoding.")
        with open(marker_genes_file, mode="r", encoding="cp1252") as file:
            reader = csv.reader(file)
            rows = list(reader)

    # Process the successfully loaded rows
    for row in rows:
        if row:
            cell_type = row[0]
            genes = row[1:]
            valid_genes = [g.strip() for g in genes if g.strip()]
            if valid_genes:
                marker_genes_dict[cell_type] = valid_genes

    available_markers = {}
    flat_gene_list = []
    seen_genes = set()  # FIXED: Eliminates overlapping markers across categories

    for cell_type, genes in marker_genes_dict.items():
        available = []
        for g in genes:
            if g in adata.var_names and g not in seen_genes:
                available.append(g)
                seen_genes.add(g)
        if available:
            available_markers[cell_type] = available
            for g in available:
                if g not in flat_gene_list:
                    flat_gene_list.append(g)

    if available_markers:
        custom_width = max(6, len(flat_gene_list) * 0.45 + 3)
        custom_height = max(4, num_clusters * 0.4 + 2)

        sc.pl.dotplot(adata, available_markers, "leiden", dendrogram=True, show=False, swap_axes=False, figsize=(custom_width, custom_height))
        plt.savefig(f'{sample_id}_marker_selectedmarkers_dotplot.pdf', bbox_inches='tight')
        plt.close()

        sc.pl.stacked_violin(adata, available_markers, groupby="leiden", swap_axes=False, dendrogram=True, show=False, figsize=(custom_width, custom_height))
        plt.savefig(f'{sample_id}_marker_selectedmarkers_stacked_violin.pdf', bbox_inches='tight')
        plt.close()

        dynamic_size = 5000 / np.sqrt(adata.n_obs)
        combined_features = ['leiden'] + flat_gene_list
        sc.pl.umap(adata, color=combined_features, s=dynamic_size, frameon=False, ncols=4, vmin=0, vmax="p99", show=False)
        plt.savefig(f'{sample_id}_marker_custom_umap_expression.pdf', bbox_inches='tight')
        plt.close()

        sc.pl.tsne(adata, color=combined_features, s=dynamic_size, frameon=False, ncols=4, vmin=0, vmax="p99", show=False)
        plt.savefig(f'{sample_id}_marker_custom_tsne_expression.pdf', bbox_inches='tight')
        plt.close()

        sc.pl.pca(adata, color=combined_features, s=dynamic_size, frameon=False, ncols=4, wspace=0.3, vmin=0, vmax="p99", sort_order=True, show=False)
        plt.savefig(f'{sample_id}_marker_custom_pca_expression.pdf', bbox_inches='tight')
        plt.close()

        sc.pl.dotplot(adata, available_markers, "leiden", dendrogram=True, show=False, swap_axes=True, figsize=(custom_width, custom_height), standard_scale='var')
        plt.savefig(f'{sample_id}_marker_selectedmarkers_dotplot_scaled.pdf', bbox_inches='tight')
        plt.close()

        

        # 1. Flatten the dictionary and verify gene presence (same as before)
        flat_marker_list = []
        for genes in available_markers.values():
            flat_marker_list.extend(genes)
        flat_marker_list = [g for g in set(flat_marker_list) if g in adata.var_names]

        # 2. DYNAMICALLY ADJUST COLUMNS BASED ON TOTAL GENES
        num_genes = len(flat_marker_list)
        if num_genes <= 2:
            ncols = num_genes
        elif num_genes <= 6:
            ncols = 3
        else:
            ncols = 4  # Max out at 4 columns for large sets of genes

        
        nrows = math.ceil(num_genes / ncols)

        # 3. DYNAMICALLY SIZE THE FIGSIZE (Give each panel a wider 5.5 inch width)
        panel_width = 5.5   # Increased from 4 to 5.5 to make panels wider/less narrow
        panel_height = 4.0  # Slightly increased height for proportion
        
        fig, axes = plt.subplots(
            nrows=nrows, 
            ncols=ncols, 
            figsize=(ncols * panel_width, nrows * panel_height)
        )
        
        # Ensure axes is a flat array even if it's a 1-panel or 1-row plot
        if num_genes == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        # 4. Loop through each gene and plot it in its own subplot panel
        for i, gene in enumerate(flat_marker_list):
            sc.pl.violin(
                adata, 
                keys=gene, 
                groupby="leiden", 
                stripplot=False, 
                ax=axes[i],      # <--- Direct the plot to the correct subplot panel
                show=False
            )
            axes[i].set_title(gene)  # Ensure the gene name is cleanly labeled on top

        # 5. Hide any blank panels at the end if your gene count doesn't perfectly fill the grid
        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        # 6. Adjust spacing and save the figure
        plt.tight_layout()
        plt.savefig(f'{sample_id}_all_markers_individual_violins.pdf', bbox_inches='tight')
        plt.close()
        

    # =========================================================================
    # 4. CELL TYPE / CLUSTER DISTRIBUTION (STACKED BAR PLOTS)
    # =========================================================================
    if target_obs:
        print("Generating stacked bar plots for cell type distributions...")
        counts = adata.obs[target_obs].value_counts()
        proportions = (counts / counts.sum() * 100).to_frame().T
        
        num_cell_types = proportions.shape[1]
        global_colors = customcolor[:num_cell_types]
        
        fig, ax = plt.subplots(figsize=(5, 8))
        proportions.plot(kind='bar', stacked=True, ax=ax, edgecolor='none', color=global_colors)
        ax.grid(False)
        ax.set_ylabel('Percentage of Sample (%)')
        ax.set_xlabel(sample_id)
        ax.set_xticklabels([]) 
        
        plt.legend(title='Cell Types', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig(f'{sample_id}_global_celltype_composition.pdf', bbox_inches='tight')
        plt.close()

        cluster_in_ct = pd.crosstab(adata.obs[target_obs], adata.obs['leiden'], normalize='index') * 100
        num_clusters_present = cluster_in_ct.shape[1]
        cluster_colors = customcolor[:num_clusters_present]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        cluster_in_ct.plot(kind='bar', stacked=True, ax=ax, edgecolor='none', color=cluster_colors)
        ax.grid(False)