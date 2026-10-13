process COMBINED_QC {
    tag "Generating Combined Cohort QC Report"
    publishDir "${params.outdir}/07_visualization", mode: 'copy'

    input:
    // Accepts raw/pre-QC data files and filtered/post-QC data files separately
    tuple val(raw_metas), path(raw_adatas)
    tuple val(qc_metas), path(qc_adatas)

    output:
    path "cohort_pre_vs_post_qc_metrics.pdf", emit: plot

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Helper function to load and format metadata dataframes
    def load_cohort_obs(files_str, metas_str, status_label):
        files = files_str.split()
        ids = metas_str.replace('[','').replace(']','').split(', ')
        
        dfs = []
        for filepath, sample_name in zip(files, ids):
            try:
                adata = sc.read_h5ad(filepath)
                df = adata.obs.copy()
                df['sample_origin'] = sample_name
                df['qc_status'] = status_label
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {filepath}: {e}")
        return dfs

    # 1. Parse and collect all dataframes
    raw_dfs = load_cohort_obs("${raw_adatas}", "${raw_metas.collect{it.id}}", "Before QC")
    qc_dfs = load_cohort_obs("${qc_adatas}", "${qc_metas.collect{it.id}}", "After QC")
    
    all_dfs = raw_dfs + qc_dfs
    if not all_dfs:
        print("Error: No data successfully parsed.")
        sys.exit(1)

    cohort_obs = pd.concat(all_dfs, axis=0)

    # Convert standard single-cell metrics to numeric type to prevent plotting failures
    metrics = ["total_counts", "pct_counts_mt", "n_genes_by_counts"]
    valid_metrics = [m for m in metrics if m in cohort_obs.columns]
    for m in valid_metrics:
        cohort_obs[m] = pd.to_numeric(cohort_obs[m], errors='coerce')

    # =========================================================================
    # GENERATE PRE- VS POST- QC COMPARISON GRID
    # =========================================================================
    if valid_metrics:
        fig, axes = plt.subplots(len(valid_metrics), 1, figsize=(12, 5 * len(valid_metrics)), sharex=False)
        if len(valid_metrics) == 1:
            axes = [axes]

        # Use Seaborn's split violin design for clean side-by-side cohort comparisons
        for i, metric in enumerate(valid_metrics):
            sns.violinplot(
                data=cohort_obs,
                x='sample_origin',
                y=metric,
                hue='qc_status',
                split=True,       # Draws left-side raw vs right-side filtered cell spreads
                inner='quart',    # Shows inner data spread markers
                palette={'Before QC': '#8A969E', 'After QC': '#0F3A85'}, # Lilly color accents
                ax=axes[i]
            )
            axes[i].set_title(f'Cohort Filter Impact - {metric}', fontsize=12, fontweight='bold')
            axes[i].set_xlabel('Sample Name')
            axes[i].set_ylabel(metric)
            axes[i].tick_params(axis='x', rotation=30)
            axes[i].grid(False)
            
            # Position the status legend cleanly outside the core plot framework
            axes[i].legend(title='Filtering Track', bbox_to_anchor=(1.02, 1), loc='upper left')

        plt.savefig('cohort_pre_vs_post_qc_metrics.pdf', bbox_inches='tight')
        plt.close()
        print("Multi-sample Pre vs Post QC distribution comparison compiled successfully.")
    else:
        print("Error: Target metric layers missing from tracking matrices.")
    """
}
