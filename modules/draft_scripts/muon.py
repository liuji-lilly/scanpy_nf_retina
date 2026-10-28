import muon as mu
import snapatac2 as snap
import scanpy as sc
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import seaborn as sns

from muon import atac as ac
from muon import prot as pt
import kaleido
import magic

# Run the built-in chrome download script
kaleido.get_chrome_sync()

# 1. Change your active working directory
os.chdir("/fsx/home/l128405/10xgenomics/multiome/S41A/")
# Verify you are in the correct place
print("Current Working Directory:", os.getcwd())

input="/fsx/home/l128405/10xgenomics/multiome/S41A/outs/filtered_feature_bc_matrix"

adata_multimodal = sc.read_10x_mtx(input)
adata_rna = adata_multimodal[:, adata_multimodal.var["feature_types"] == "Gene Expression"].copy()
adata_rna.var_names_make_unique()

# Save to h5ad
adata_rna.write("/fsx/home/l128405/10xgenomics/multiome/S41A/filtered_feature_bc_matrix_rna_only.h5ad")

# Read via the MTX Folder Layout
#mdata = mu.read_10x_mtx("/fsx/home/l128405/10xgenomics/multiome/S41A/outs/filtered_feature_bc_matrix")
#print(mdata)
# 1. Initialize the backed AnnData object from the fragment file
# Use the mm39 genome assembly for fragment chromosome parsing
# Define your input paths clearly
fragment_path = "/fsx/home/l128405/10xgenomics/multiome/S41A/outs/atac_fragments.tsv.gz"
output_file = "/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_preprocessed.h5ad"
# Ensure paths align with your active directory setup
pdf_output_path = "/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_qc_report.pdf"

genome_assembly = snap.genome.mm39 

# Ensure old files are purged before execution
if os.path.exists(output_file):
    os.remove(output_file)

data_atac = snap.pp.import_fragments(
    fragment_file=fragment_path,
    chrom_sizes=genome_assembly, 
    file=output_file, 
    is_paired=True,
    sorted_by_barcode=False,
    min_num_fragments=400
)

print(data_atac)


# 2. Compute TSS enrichment scores
snap.metrics.tsse(data_atac, gene_anno=genome_assembly)

# 3. Plot metrics if needed, or filter directly using strict thresholds
# Initialize the multi-page PDF stream context
# Define your genome annotation (e.g., snap.genome.mm39 or hg38 depending on your data)
# Make sure to replace this with the genome assembly you are actively using

with PdfPages(pdf_output_path) as pdf:
    # ----------------------------------------------------
    # Page 1: Standard Fragment Size Distribution (Linear)
    # ----------------------------------------------------
    print("Computing and writing Page 1: Linear Fragment Size Distribution...")
    # Calculate the raw vector directly (returns a 1D numpy array)
    frag_dist = snap.metrics.frag_size_distr(data_atac, inplace=False)
    
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(frag_dist, color="darkcyan", lw=1.5)
    ax.set_title("Fragment Size Distribution")
    ax.set_xlabel("Fragment Size (bp)")
    ax.set_ylabel("Count")
    ax.set_xlim(0, 1000)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # ----------------------------------------------------
    # Page 2: Log-Y Fragment Size Distribution
    # ----------------------------------------------------
    print("Writing Page 2: Log-Y Fragment Size Distribution...")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(frag_dist, color="darkcyan", lw=1.5)
    ax.set_yscale("log")
    ax.set_title("Fragment Size Distribution (Log Scale)")
    ax.set_xlabel("Fragment Size (bp)")
    ax.set_ylabel("Count (Log)")
    ax.set_xlim(0, 1000)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # ----------------------------------------------------
    # Page 3: TSS Enrichment Score Panel
    # ----------------------------------------------------
    print("Computing and writing Page 3: TSS Enrichment Evaluation Plot...")
    # Compute the raw TSSE profiles dictionary
    # Calculate cell-level TSSE metrics and save them directly into data_atac.obs
    snap.metrics.tsse(data_atac, genome_assembly, inplace=True)
    
    # Extract the pre-calculated metrics from the metadata fields
    n_frag = data_atac.obs['n_fragment']
    tsse_score = data_atac.obs['tsse']
    
    # Initialize the matplotlib figure canvas
    fig, ax = plt.subplots(figsize=(7, 5.5))
    
    # Plot the 2D density contour map matching the Plotly architecture
    sns.kdeplot(
        x=n_frag,
        y=tsse_score,
        fill=True,
        thresh=0.002,          # Sets the outer boundary threshold matching the colorbar
        levels=10,             # Number of concentric contour bands
        cmap="Blues",          # Matching the blue color palette layout
        cbar=True,             # Display the vertical color scale bar
        cbar_kws={'label': 'Density'},
        log_scale=(True, False), # Strictly logs the X-axis (fragments) while leaving Y linear
        ax=ax
    )
    
    # Adjust titles and layout properties to fit the standard single-cell reporting guidelines
    ax.set_title("TSS Enrichment vs. Unique Fragments QC", fontsize=12, pad=10)
    ax.set_xlabel("Number of unique fragments", fontsize=10)
    ax.set_ylabel("TSS enrichment score", fontsize=10)
    
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)    



# Standard 10x multiome thresholds: TSSE > 5, fragments >2000
snap.pp.filter_cells(data_atac, min_counts=2000, max_counts=100000, min_tsse=5)


# 4. Generate a genome-wide 500bp tile matrix
snap.pp.add_tile_matrix(data_atac, bin_size=500)

# 5. Select highly variable features across chromatin tiles
snap.pp.select_features(data_atac, n_features=150000)

# 6. Doublet removal
snap.pp.scrublet(data_atac)
snap.pp.filter_doublets(data_atac)

output_file = "/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_qc.h5ad"
data_atac.write(output_file)

genome_assembly = snap.genome.mm39 

import pandas as pd

# 7. Dimension reduction 
snap.tl.spectral(data_atac)

# 8. Clustering analysis: Build k-nearest neighbor graph using the embedding vectors
snap.pp.knn(data_atac, use_rep="X_spectral", n_neighbors=15)
snap.tl.leiden(data_atac, key_added="atac_leiden", resolution=2)
snap.tl.umap(data_atac)

snap.pl.umap(data_atac, color=["atac_leiden", "n_fragment"], show=False, out_file="S41A_atac_umap.html", height=500)
snap.pl.umap(data_atac, color=["atac_leiden", "n_fragment"], show=False, out_file="S41A_atac_umap.pdf", height=500)

# Create and sort the distribution dataframe
cluster_df = data_atac.obs['atac_leiden'].value_counts().to_pandas()
cluster_df['atac_leiden'] = cluster_df['atac_leiden'].astype(int)
cluster_df = cluster_df.sort_values("atac_leiden").reset_index(drop=True)

# Add percentage tracking for your workflow metrics
total_cells = cluster_df['count'].sum()
cluster_df['percentage (%)'] = (cluster_df['count'] / total_cells * 100).round(2)

cluster_df.to_csv('S41A_cluster_sizes.csv', index=False)

 
output_file = "/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_clustered.h5ad"
data_atac.write(output_file)


#  Imputation: MAGIC algorithm to perform imputation and data smoothing
input_file = "/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_clustered.h5ad"
marker_genes = ['Myo7a', 'Sox2', 'Tmc1', 'Tmc2']

genome_assembly = snap.genome.mm39 

data_atac = snap.read(input_file)
gene_matrix = snap.pp.make_gene_matrix(data_atac, genome_assembly)
sc.pp.filter_genes(gene_matrix, min_cells= 5)
sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)
sc.external.pp.magic(gene_matrix, solver="approximate")
gene_matrix.obsm["X_umap"] = data_atac.obsm["X_umap"]

# Generate the plot structure without rendering it to screen
sc.pl.umap(
    gene_matrix, 
    use_raw=False, 
    color=["atac_leiden", "n_fragment"] + marker_genes,
    show=False
)
# Export the active figure to any absolute path on your cluster
plt.savefig(
    "/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_umap_marker_genes.pdf", 
    bbox_inches='tight', 
    dpi=150
)
plt.close()

gene_matrix.write("/fsx/home/l128405/10xgenomics/multiome/S41A/S41A_atac_mat.h5ad", compression='gzip')

flat_gene_list=marker_genes
num_clusters=len(data_atac.obs['atac_leiden'].unique())
custom_width = max(6, len(flat_gene_list) * 0.45 + 3)
custom_height = max(4, num_clusters * 0.4 + 2)
custom_width_violin = max(15, len(flat_gene_list) * 0.6 + 4)
custom_height_violin = max(18, num_clusters * 0.6 + 3) # Significantly increased from 0.4

sc.pl.dotplot(
    gene_matrix, 
    flat_gene_list, 
    groupby="atac_leiden",  # Use your high-res clustering key
    dendrogram=True, 
    show=False, 
    swap_axes=False, 
    figsize=(custom_width, custom_height)
)

plt.savefig('S41A_marker_selectedmarkers_dotplot.pdf', bbox_inches='tight')
plt.close()

sc.pl.stacked_violin(gene_matrix, flat_gene_list, groupby="atac_leiden", swap_axes=False, dendrogram=True, show=False,dodge=False, figsize=(custom_width_violin, custom_height_violin))
plt.savefig('S41A_marker_selectedmarkers_stacked_violin.pdf', bbox_inches='tight')
plt.close()


# 9. Call cluster-specific peaks safely
snap.tl.macs3(data_atac, groupby="atac_leiden")

# 10. Merge single cluster peak coordinates into a consensus peak matrix
peaks = snap.tl.merge_peaks(data_atac.uns['macs3'], genome_assembly)

# 11. Project raw fragment counts into the new clean consensus peak matrix
peak_mat = snap.pp.make_peak_matrix(data_atac, use_rep=peaks['Peaks'])

#12. Finding marker regions
marker_peaks = snap.tl.marker_regions(peak_mat, groupby='atac_leiden', pvalue=0.01)
snap.pl.regions(peak_mat, groupby='atac_leiden', peaks=marker_peaks, interactive=False,out_file="S41A_actac_leiden_marker_regions.pdf", height=500)

# 12. Finding marker regions and generating the plot structure
fig = snap.pl.regions(
    peak_mat, 
    groupby='atac_leiden', 
    peaks=marker_peaks, 
    interactive=False, 
    height=600  # Give the canvas a bit more breathing room
)

# Force Plotly to display every single cluster label on the X-axis
fig.update_xaxes(
    tickmode='linear',
    dtick=1,
    tickangle=45  # Rotate the text labels to prevent side-by-side collisions
)

# Export the clean, fully-labeled figure to your target path
fig.write_image("S41A_actac_leiden_marker_regions.pdf")


motifs = snap.tl.motif_enrichment(
    motifs=snap.datasets.cis_bp(unique=True),
    regions=marker_peaks,
    genome_fasta=genome_assembly,
)
snap.pl.motif_enrichment(motifs, max_fdr=0.0001, height=1600, interactive=False,out_file="S41A_actac_leiden_motif_enrichment.pdf")
import pandas as pd

# Combine everything cleanly into one Excel workbook
with pd.ExcelWriter("S41A_actac_leiden_motif_enrichment_results.xlsx") as writer:
    for cluster_name, df in motifs.items():
        # df.to_pandas() is natively optimized in Polars
        df.to_pandas().to_excel(writer, sheet_name=f"Cluster_{cluster_name}", index=False)

for cluster_name, df in motifs.items():
    df.write_csv(f"S41A_actac_motif_enrichment_cluster_{cluster_name}.csv")

import polars as pl
# 1. Identify the unique top motifs globally across all clusters
top_motif_ids = set()
for cluster, df in motifs.items():
    sig_df = df.filter(pl.col("adjusted p-value") < 0.05).sort("log2(fold change)", descending=True)
    cluster_top_10 = sig_df.head(10)["id"].to_list()
    top_motif_ids.update(cluster_top_10)

# 2. Dynamic Height Calculation
num_motifs = len(top_motif_ids)
dynamic_height = max(600, (num_motifs * 14) + 300) 
print(f"Total unique motifs to plot: {num_motifs}. Setting plot height to: {dynamic_height}px")

# 3. Filter the original data dict and shorten the text to prevent overlapping
filtered_motifs = {}
for cluster, df in motifs.items():
    filtered_df = (
        df.filter(pl.col("id").is_in(list(top_motif_ids)))
          .sort("id")
          .with_columns(pl.col("id").str.split("+").list.get(0)) 
    )
    filtered_motifs[cluster] = filtered_df

# 4. Generate the Figure object directly (remove out_file here)
fig = snap.pl.motif_enrichment(
    filtered_motifs,
    show=False,
    height=dynamic_height,  
    width=1100             
)

# 5. FORCE CRISP PIXELS: Disable smooth blending on the heatmap trace
fig.update_traces(
    zsmooth=False,                  # Disables Plotly's internal blurring/interpolation
    selector=dict(type='heatmap')
)
fig.write_html("S41A_actac_leiden_motif_enrichment_top10.html")

# 6. Save the crisp canvas using an explicit scale multiplier
# scale=3 increases the base resolution layout to force sharp grid lines in PDF viewers
fig.write_image(
    "S41A_actac_leiden_motif_enrichment_top10.pdf", 
    engine="kaleido", 
    scale=3
)

# 12. Read your matched gene expression matrix 
# (Parsing MTX layout avoids the known Unicode h5py decoding crashes)
data_rna = sc.read_10x_mtx("/fsx/home/l128405/10xgenomics/multiome/S41A/S41A/outs/filtered_feature_bc_matrix")

# Ensure cell barcodes match explicitly across both structures
shared_barcodes = list(set(data_rna.obs_names) & set(data_peaks.obs_names))
data_rna = data_rna[shared_barcodes].copy()
data_peaks = data_peaks[shared_barcodes].copy()

# 13. Create your unified MuData object
mdata = mu.MuData({
    "rna": data_rna,
    "atac": data_peaks
})

print("🚀 Joint Multiome Object Built Successfully:")
print(mdata)


