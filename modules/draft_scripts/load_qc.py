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

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white')

# Load data
print("Loading 10X data...")
adata = sc.read_10x_mtx('${input_dir}', var_names='gene_symbols', cache=True)

#adata = sc.read_10x_mtx('/fsx/home/l128405/scanpy_nf_retina/sampledata/test1', var_names='gene_symbols', cache=True)

# Initial stats
n_cells_initial = adata.n_obs
n_genes_initial = adata.n_vars

# mitocondrial genes
adata.var['mt'] =adata.var_names.str.startswith('mt-')| adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('MT-')
# ribosomal genes
#adata.var["ribo"] = adata.var_names.str.startswith(("rps", "rpl"))  | adata.var_names.str.startswith(("Rps", "Rpl"))  | adata.var_names.str.startswith(("RPS", "RPL")) 
# hemoglobin genes.
#adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")| adata.var_names.str.contains(r"^Hb[abdegmqz]\d*(?!\w)")| adata.var_names.str.contains(r"^hb[abdegmqz]\d*(?!\w)")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, percent_top=[20], log1p=True
)

median_genes_initial = adata.obs['n_genes_by_counts'].median()


# --- GENERATE PRE-FILTER FIGURE ---
# Define a wide canvas for the three plots
fig_pre, axes_pre = plt.subplots(1, 3, figsize=(15, 5))

# Loop through and plot each metric onto its designated subplot axis
metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
for ax, metric in zip(axes_pre, metrics):
    sc.pl.violin(
        adata, 
        keys=metric, 
        jitter=0.4, 
        ax=ax, 
        show=False, 
        inner='quartiles'  
    )

fig_pre.suptitle(f"Pre-Filtering QC (Cells: {n_cells_initial})", y=1.02, fontsize=14)
plt.tight_layout()

# Force a clean canvas for subsequent filtering and plotting
plt.figure()



# Filter
print("Filtering cells...")
sc.pp.filter_genes(adata, min_cells=${min_cells})

#sc.pp.filter_genes(adata, min_cells=3)
#sc.pp.filter_cells(adata, min_genes=200)

# Mark outliers based on QC metrics using Median Absolute Deviation (MAD)
#QC covariates each with a threshold of 5 MADs
# --- STEP 0: DEFINE CUSTOM BIOLOGY-DRIVEN OVERRIDES ---
# (Uses your existing max_genes and max_mt_percent variables)

caps_and_floors = {
    "log1p_total_counts":        {"min": np.log1p(${min_transcripts}),   "max": np.log1p(${max_transcripts})},
    "log1p_n_genes_by_counts":   {"min": np.log1p(${min_genes}),   "max": np.log1p(${max_genes})},
    "pct_counts_in_top_20_genes": {"min": 0.0,           "max": ${max_pct_counts_in_top_20_genes}},
    "pct_counts_mt":              {"min": 0.0,            "max": ${max_mt_percent}}
}

# caps_and_floors = {
#     "log1p_total_counts":          {"min": np.log1p(500),  "max": np.log1p(15000)},  # Log thresholds for 500 & 10,000 transcripts
#     "log1p_n_genes_by_counts":     {"min": np.log1p(200),  "max": np.log1p(8000)},   # Log thresholds for 200 & 3,000 genes
#     "pct_counts_in_top_20_genes":  {"min": 0.0,            "max": 80.0},
#     "pct_counts_mt":               {"min": 0.0,            "max": 5.0}
# }

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

#adata.obs.outlier.value_counts()
# --- STEP 2: HANDLE MITOCHONDRIAL METRIC DYNAMICALLY ---
mt_data = adata.obs["pct_counts_mt"].to_numpy()
mt_median = np.median(mt_data)
mt_mad = median_abs_deviation(mt_data, scale=1.0)
mt_mad = 1e-6 if mt_mad == 0 else mt_mad

# Apply 3-MAD fence for MT, capped at max_mt_percent
mt_upper_fence = mt_median + (3 * mt_mad)
final_mt_upper = np.maximum(mt_upper_fence, caps_and_floors["pct_counts_mt"]["max"])

# Save mitochondrial outlier flag back to AnnData object
adata.obs["mt_outlier"] = mt_data > final_mt_upper

qc_summary_data.append({
    "Metric": "pct_counts_mt",
    "Median": mt_median,
    "MAD": mt_mad,
    "Lower_Fence": np.nan,
    "Upper_Fence": mt_upper_fence,
    "Final_Lower": caps_and_floors["pct_counts_mt"]["min"],
    "Final_Upper": final_mt_upper
})
##adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 20)
#adata.obs.mt_outlier.value_counts()
pd.DataFrame(qc_summary_data).round(1).to_csv('qc_cutoffs.csv', index=False)


print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

n_cells_highquality = adata.n_obs
median_genes_highquality = adata.obs['n_genes_by_counts'].median()
n_genes_highquality = adata.n_vars

print("Running scDblFinder doublet detection...")

# Run scDblFinder directly as a single standalone step
sdf = ScDblFinder(adata, random_state=0)
#sdf.run(dbr=0.03)                            # 7% expected doublet rate
sdf.run(dbr=${expected_doublet_rate})                            # 7% expected doublet rate
# Map the outputs directly to your single-cell metadata table
adata.obs.scDblFinder_class.value_counts()
# Count how many doublets were found
n_doublets = (adata.obs["scDblFinder_class"] == "doublet").sum()
n_singlets = (adata.obs["scDblFinder_class"] == "singlet").sum()
median_genes_singlets = adata.obs[adata.obs["scDblFinder_class"] == "singlet"]['n_genes_by_counts'].median()
n_genes_singlets = adata[adata.obs["scDblFinder_class"] == "singlet"].n_vars

print(f"Detected {n_doublets} doublets out of {adata.n_obs} cells.")

data = {
    "Cell Count": [n_cells_initial, n_cells_highquality, n_singlets],
    "Gene Count": [n_genes_initial, n_genes_highquality, n_genes_singlets],
    "Median Genes per Cell": [median_genes_initial, median_genes_highquality, median_genes_singlets]
}

# Save QC metrics
stages = ["Initial", "QC", "Singlet"]
qc_df = pd.DataFrame.from_dict(data, orient='index', columns=stages)
qc_df.to_csv('qc_metrics.csv', index=False)

# QC violin plot
# --- STEP 2: Generate and save Post-filter plot to a variable ---
fig_post, axes_post = plt.subplots(1, 3, figsize=(15, 5))

for ax, metric in zip(axes_post, metrics):
    sc.pl.violin(
        adata,  # Using your filtered/active object
        keys=metric, 
        jitter=0.4, 
        ax=ax, 
        show=False,inner='quartiles'  
    )

fig_post.suptitle(f"Post-Filtering QC (Cells: {adata.n_obs})", y=1.02, fontsize=14)
fig_post.tight_layout()


# --- STEP 3: Compile them into the PDF at the end of the script ---
with PdfPages('qc_plots.pdf') as pdf:
    pdf.savefig(fig_pre, bbox_inches='tight')
    pdf.savefig(fig_post, bbox_inches='tight')

# Clean up all figure memory
plt.close(fig_pre)
plt.close(fig_post) 
    
# Save
adata.write('adata_qc.h5ad')
print(f"QC complete: {adata.n_obs} cells, {adata.n_vars} genes")