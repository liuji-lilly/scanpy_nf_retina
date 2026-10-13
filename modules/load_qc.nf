process LOAD_AND_QC {
    tag "QC filtering: ${meta.id}"
    publishDir "${params.outdir}/01_qc", mode: 'copy'
    
    input:
    tuple val(meta), path(input_dir_or_file) // Only take the data channel—remove all val() items below this!

    output:
    tuple val(meta), path("${meta.id}_raw.h5ad"),        emit: adata_raw
    tuple val(meta), path("${meta.id}_qc.h5ad"),         emit: adata_qc
    tuple val(meta), path("${meta.id}_qc_cutoffs.csv"),  emit: qc_cutoff_summary
    tuple val(meta), path("${meta.id}_qc_metrics.csv"),  emit: qc_metrics
    tuple val(meta), path("${meta.id}_qc_plots.pdf"),    emit: qc_plot

    script:
    def chem = meta.chemistry
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import median_abs_deviation
    from scipy.sparse import csc_matrix
    import omicverse as ov
    from pyscdblfinder import ScDblFinder
    from matplotlib.backends.backend_pdf import PdfPages
    import os
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white')

    # Load data
    input_path = '${input_dir_or_file}'
    print(f"Loading data for ${meta.id} from: {input_path}")
    
    if os.path.isdir(input_path):
        print("Detected directory structure. Reading as Cell Ranger MTX format...")
        adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', cache=True)
    elif input_path.endswith('.h5ad'):
        print("Detected AnnData H5AD file. Reading via sc.read_h5ad...")
        adata = sc.read_h5ad(input_path)
        # If the file came from an old cellbender run and has ensembl IDs in index, resolve it:
        if 'gene_name' in adata.var.columns:
            adata.var['ensembl_id'] = adata.var_names
            adata.var_names = adata.var['gene_name']
        adata.var_names_make_unique()
    elif input_path.endswith('.h5') or input_path.endswith('.hdf5'):
        print("Detected H5 file. Using CellBender loader to restore MGI Gene Symbols...")
        from cellbender.remove_background.downstream import anndata_from_h5
        
        # FIXED: This automatically sets var_names to Gene Symbols and keeps counts structured
        adata = anndata_from_h5(input_path)
        adata.var_names_make_unique()
    else:
        raise ValueError(f"Unsupported file format found for input: {input_path}")
        
    # INSERT THE ROUNDING PROTOCOL HERE
    if 'cellbender_filter' in input_path:
        print("Rounding non-integer CellBender background probabilities back to whole counts...")
        from scipy.sparse import issparse
        if issparse(adata.X):
            adata.X.data = np.round(adata.X.data)
        else:
            adata.X = np.round(adata.X)
    
    # Initial stats
    n_cells_initial = adata.n_obs
    n_genes_initial = adata.n_vars

    # Mitochondrial genes identification
    adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('MT-')

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, percent_top=[20], log1p=True
    )

    median_genes_initial = adata.obs['n_genes_by_counts'].median()

    # --- GENERATE PRE-FILTER FIGURE ---
    fig_pre, axes_pre = plt.subplots(1, 3, figsize=(15, 5))
    metrics_list = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
    for ax, metric in zip(axes_pre, metrics_list):
        sc.pl.violin(adata, keys=metric, jitter=0.4, ax=ax, show=False, inner='quartiles')

    fig_pre.suptitle(f"Pre-Filtering QC (Cells: {n_cells_initial})", y=1.02, fontsize=14)
    plt.tight_layout()
    plt.figure()

    # Filter out empty genes natively
    print("Filtering cells...")
    sc.pp.filter_genes(adata, min_cells=${params.min_cells ?: 3})

    # Global params injection mapping with structural fallbacks
    caps_and_floors = {
        "log1p_total_counts":        {"min": np.log1p(${params.min_transcripts ?: 500}),"max": np.log1p(${params.max_transcripts ?: 10000})},
        "log1p_n_genes_by_counts":   {"min": np.log1p(${params.min_genes ?: 200}),"max": np.log1p(${params.max_genes ?: 3000})},
        "pct_counts_in_top_20_genes": {"min": 0.0,"max": ${params.max_pct_counts_in_top_20_genes ?: 80.0}},
        "pct_counts_mt":              {"min": 0.0,"max": ${params.max_mt_percent ?: 5.0}}
    }

    # --- STEP 1: COMPUTE BOUNDARIES FOR MAIN QC MATRIX ---
    metrics = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes"]
    data_matrix = adata.obs[metrics].to_numpy()

    # Compute medians and MADs simultaneously across all columns
    medians = np.median(data_matrix, axis=0)
    mads = median_abs_deviation(data_matrix, axis=0, scale=1.0)

    # Avoid divide-by-zero errors on ultra-uniform columns
    mads = np.where(mads == 0, 1e-6, mads)

    # Calculate raw mathematical 5-MAD fences
    lower_fences = medians - (5 * mads)
    upper_fences = medians + (5 * mads)

    # Dynamic clamping loop to tag outliers
    mad_outliers = np.zeros_like(data_matrix, dtype=bool)
    qc_summary_data = []
    for i, metric in enumerate(metrics):
        limits = caps_and_floors[metric]
        
        # Use the caps and floor values directly as the final filtering cutoffs
        final_lower = limits["min"]
        final_upper = limits["max"]
        
        # Flag the outliers based on the hard final cutoffs
        mad_outliers[:, i] = (data_matrix[:, i] < final_lower) | (data_matrix[:, i] > final_upper)
        
        # Conditionally check if the metric is log-transformed
        is_log = metric.startswith("log1p_")
        if is_log:
            # Convert log1p values back into raw count integers
            med_raw     = np.expm1(medians[i])
            low_raw     = max(0, np.expm1(lower_fences[i]))
            mad_raw     = np.expm1(medians[i] + mads[i]) - med_raw 
            upp_raw     = np.expm1(upper_fences[i])
            f_low_raw   = max(0, np.expm1(final_lower))
            f_upp_raw   = np.expm1(final_upper)
        else:
            # Keep percentages (e.g. pct_counts_in_top_20_genes) as-is
            med_raw     = medians[i]
            low_raw     = lower_fences[i]
            mad_raw     = mads[i]
            upp_raw     = upper_fences[i]
            f_low_raw   = final_lower
            f_upp_raw   = final_upper
        
        qc_summary_data.append({
            "Metric": metric,
            "Median": med_raw,
            "MAD": mad_raw,
            "Lower_Fence": low_raw,
            "Upper_Fence": upp_raw,
            "Final_Lower": f_low_raw,
            "Final_Upper": f_upp_raw
        })


    # Save main outlier flag back to AnnData object
    adata.obs["outlier"] = np.any(mad_outliers, axis=1)

    # --- STEP 2: HANDLE MITOCHONDRIAL METRIC DYNAMICALLY ---
    mt_data = adata.obs["pct_counts_mt"].to_numpy()
    mt_median = np.median(mt_data)
    mt_mad = median_abs_deviation(mt_data, scale=1.0)
    mt_mad = 1e-6 if mt_mad == 0 else mt_mad

    mt_upper_fence = mt_median + (3 * mt_mad)
    final_mt_upper = np.maximum(mt_upper_fence, caps_and_floors["pct_counts_mt"]["max"])
    adata.obs["mt_outlier"] = mt_data > final_mt_upper

    qc_summary_data.append({
        "Metric": "pct_counts_mt", "Median": mt_median, "MAD": mt_mad,
        "Lower_Fence": np.nan, "Upper_Fence": mt_upper_fence,
        "Final_Lower": caps_and_floors["pct_counts_mt"]["min"], "Final_Upper": final_mt_upper
    })
    pd.DataFrame(qc_summary_data).round(1).to_csv('${meta.id}_qc_cutoffs.csv', index=False)

    # Save checkpoint raw data asset before subset filter
    adata.write('${meta.id}_raw.h5ad')

    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    n_cells_highquality = adata.n_obs
    median_genes_highquality = adata.obs['n_genes_by_counts'].median()
    n_genes_highquality = adata.n_vars

    print("Running scDblFinder doublet detection...")
    
    # Dynamic calculation matching the loaded single-cell chemistry choice
    chemistry_type = "${chem}".upper()
    if "GEM-X" in chemistry_type:
        # GEM-X rate: ~0.4% per 1,000 cells recovered
        target_dbr = (adata.n_obs / 1000.0) * 0.004
    else:
        # Standard Next GEM rate: ~0.8% per 1,000 cells recovered
        target_dbr = (adata.n_obs / 1000.0) * 0.008

    print(f"Calculated expected doublet rate: {target_dbr:.4f} based on chemistry: {chemistry_type}")

    sdf = ScDblFinder(adata, random_state=0)
    sdf.run(dbr=target_dbr)                            
    
    adata.obs['predicted_doublet'] = adata.obs['scDblFinder_class'] == 'doublet'
    n_doublets = int((adata.obs["scDblFinder_class"] == "doublet").sum())
    n_singlets = int((adata.obs["scDblFinder_class"] == "singlet").sum())
    
    median_genes_singlets = adata.obs[adata.obs["scDblFinder_class"] == "singlet"]['n_genes_by_counts'].median()
    n_genes_singlets = adata[adata.obs["scDblFinder_class"] == "singlet"].n_vars

    print(f"Detected {n_doublets} doublets out of {adata.n_obs} cells.")
    
    # Subset to keep pure singlets matrix
    adata = adata[adata.obs["scDblFinder_class"] == "singlet"].copy()

    # Clean multi-stage matrix tracking with sample_id prepended as first column row
    qc_df = pd.DataFrame(
        [["${meta.id}", n_cells_initial, n_genes_initial, median_genes_initial, n_cells_highquality, n_genes_highquality, median_genes_highquality, n_singlets, n_genes_singlets, median_genes_singlets, n_doublets]], 
        columns=['sample_id', '#cells_initial', '#genes_initial', '#median_genes_initial', '#cells_post_mad', '#genes_post_mad', '#median_genes_post_mad', '#cells_final_singlets', '#genes_final_singlets', '#median_genes_final_singlets', '#doublets_removed']
    )
    qc_df.to_csv('${meta.id}_qc_metrics.csv', index=False)

    # QC violin plot generation step
    fig_post, axes_post = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric in zip(axes_post, metrics_list):
        sc.pl.violin(adata, keys=metric, jitter=0.4, ax=ax, show=False, inner='quartiles')

    fig_post.suptitle(f"Post-Filtering QC (Cells: {adata.n_obs})", y=1.02, fontsize=14)
    fig_post.tight_layout()

    # Compile images into multi-page diagnostic PDF channel asset target
    with PdfPages('${meta.id}_qc_plots.pdf') as pdf:
        pdf.savefig(fig_pre, bbox_inches='tight')
        pdf.savefig(fig_post, bbox_inches='tight')

    plt.close(fig_pre)
    plt.close(fig_post) 
        
    # Final active output save matrix
    adata.write('${meta.id}_qc.h5ad')
    print(f"QC complete for ${meta.id}: {adata.n_obs} cells, {adata.n_vars} genes remaining.")
    """
}