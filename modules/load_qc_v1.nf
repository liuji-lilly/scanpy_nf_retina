process LOAD_AND_QC {
    tag "QC filtering: ${meta.id}"
    publishDir "${params.outdir}/01_qc", mode: 'copy'
    
    input:
    tuple val(meta), path(input_dir)
    val(min_genes)
    val(min_cells)
    val(max_genes)
    val(max_mt_percent)
    
    output:
    tuple val(meta), path("${meta.id}_raw.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_qc.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_qc_metrics.csv"), emit: qc_metrics
    tuple val(meta), path("${meta.id}_qc_plots.pdf"), emit: qc_plot
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white')
    
    # Load data
    print(f"Loading 10X data for ${meta.id}...")
    adata = sc.read_10x_mtx('${input_dir}', var_names='gene_symbols', cache=True)
    
    # INJECT METADATA INTO THE ANNDATA OBJECT
    # This stores your metadata in the cell-level table (.obs)
    adata.obs['sample_id'] = "${meta.id}"
    adata.obs['batch']     = "${meta.batch}"
    adata.obs['treatment'] = "${meta.treatment}"
    
    # Initial stats
    n_cells_initial = adata.n_obs
    n_genes_initial = adata.n_vars
    
    # Calculate QC metrics to get n_genes_by_counts
    adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                                log1p=False, inplace=True)
                                
    median_genes_initial = adata.obs['n_genes_by_counts'].median()

    # --- GENERATE PRE-FILTER FIGURE ---
    metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    fig_pre, axes_pre = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric in zip(axes_pre, metrics):
        sc.pl.violin(adata, keys=metric, jitter=0.4, ax=ax, show=False)
    
    fig_pre.suptitle(f"Pre-Filtering QC (Cells: {int(n_cells_initial)})", y=1.02, fontsize=14)
    fig_pre.tight_layout()

    adata.write('${meta.id}_raw.h5ad')

    # Filter cells and genes
    print("Filtering cells...")
    sc.pp.filter_cells(adata, min_genes=${min_genes})
    sc.pp.filter_genes(adata, min_cells=${min_cells})
    
    print("Running Scrublet doublet detection...")
    sc.pp.scrublet(adata)
    
    # Count how many doublets were found
    n_doublets = sum(adata.obs['predicted_doublet'])
    print(f"Detected {n_doublets} doublets out of {adata.n_obs} cells.")
    
    # Filter out the doublets
    adata = adata[~adata.obs['predicted_doublet'], :].copy()

    # Apply feature cutoff filters
    adata = adata[adata.obs.n_genes_by_counts < ${max_genes}, :]
    adata = adata[adata.obs.pct_counts_mt < ${max_mt_percent}, :]
    
    # Save QC metrics using the dynamic sample name
    qc_df = pd.DataFrame({
        'metric': ['#cells', '#genes', '#median genes/cell', '#cells post filtering', '#unique genes post filtering','#median genes/cell post filtering', '#doublets_removed'],
        'value': [n_cells_initial, n_genes_initial, median_genes_initial, adata.n_obs, adata.n_vars, adata.obs['n_genes_by_counts'].median(), n_doublets]
    })
    qc_df.to_csv('${meta.id}_qc_metrics.csv', index=False)
    
    # --- GENERATE POST-FILTER FIGURE ---
    fig_post, axes_post = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric in zip(axes_post, metrics):
        sc.pl.violin(adata, keys=metric, jitter=0.4, ax=ax, show=False)
        
    fig_post.suptitle(f"Post-Filtering QC (Cells: {adata.n_obs})", y=1.02, fontsize=14)
    fig_post.tight_layout()

    # --- COMPILE THEM INTO THE PDF ---
    with PdfPages('${meta.id}_qc_plots.pdf') as pdf:
        pdf.savefig(fig_pre, bbox_inches='tight')
        pdf.savefig(fig_post, bbox_inches='tight')

    # Clean up figure memory
    plt.close(fig_pre)
    plt.close(fig_post) 
        
    # Save AnnData with unique filename
    adata.write('${meta.id}_qc.h5ad')
    print(f"QC complete for ${meta.id}: {adata.n_obs} cells, {adata.n_vars} genes")
    """
}
