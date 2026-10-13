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
sys.path.append(os.path.join("/fsx/home/l128405/scanpy_nf_retina", "bin"))
from sctypepy.sctype_py import gene_sets_prepare, sctype_score, process_cluster

sc.settings.verbosity = 3

# ==========================================
# 1. READ ANNDATA FIRST
# ==========================================
print("Loading dimensionality reduction data for test1...")
adata = sc.read_h5ad("/fsx/home/l128405/scanpy_nf_retina/sampledata/results/06_annotation/test1_adata_annotated.h5ad")

# ==========================================
# 2. READ CSV & EXTRACT CELLS OF INTEREST
# ==========================================
print("Reading target cells of interest from CSV...")
df_interest = pd.read_csv('/fsx/home/l128405/scanpy_nf_retina/sampledata/subclustering.csv')

# Clean header whitespaces/tabs if any exist
df_interest.columns = df_interest.columns.str.strip()

# Handle the fact that your CSV says "test2" but meta.id says "test2_adata"
clean_id = "test1".replace("_adata", "").strip()
# Locate column 0 dynamically rather than using strict strings to dodge naming glitches
id_col_name = df_interest.columns[0]
df_sample = df_interest[df_interest[id_col_name].astype(str).str.strip() == clean_id]

if df_sample.empty:
    raise ValueError(f"Error: Sample ID '{clean_id}' derived from 'test1' not found in CSV column '{id_col_name}'!")

# Extract target properties safely from rows
target_column_in_obs = df_sample.iloc[0, 1].strip()          # Extracts 'sctype_classification'
allowed_identities = [df_sample.iloc[0, 2].strip()]          # Extracts ['Hair Cells']

# Create a clean file-system-safe identifier string from the identities list
identities_str = "_".join([str(x) for x in allowed_identities])
identities_safe = re.sub(r'[^a-zA-Z0-9_\-]', '', identities_str.replace(' ', '_'))

# Prefix identifier structure for outputs
file_prefix = f"test1.{identities_safe}"
print(f"File names for this run will use prefix: {file_prefix}")

print(f"Filtering cells where '{target_column_in_obs}' is in: {allowed_identities}")

if target_column_in_obs not in adata.obs.columns:
    raise ValueError(f"Column '{target_column_in_obs}' not found in adata.obs metadata!")

# Subset the AnnData object to only keep cells matching your criteria
adata = adata[adata.obs[target_column_in_obs].isin(allowed_identities)].copy()
print(f"Subset complete. Remaining cells: {adata.n_obs}")

if adata.n_obs == 0:
    raise ValueError("Error: Zero cells left after filtering down with interest_csv!")

# ==========================================
# 3. RERUN PIPELINE (HVG -> DIM RED -> CLUSTER)
# ==========================================
print("Identify HVG...")
sc.pp.highly_variable_genes(adata, layer="scran_normalization")

plt.figure(figsize=(6, 5))
ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=15
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.savefig(f'{file_prefix}_hvg_plot.pdf', format="pdf", bbox_inches="tight", dpi=300)
plt.close()

print("Computing PCA...")
adata.var["highly_variable"] = adata.var["highly_deviant"]
sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable")

# Plot variance explained
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
fig = plt.gcf()
fig.suptitle(f"PCA Variance Ratio - test1", y=1.02, fontsize=12)
plt.savefig(f'{file_prefix}_pca_variance.pdf', format="pdf", bbox_inches="tight", dpi=300)
plt.close()

print("Computing neighborhood graph...")
# Swapped from global params to local input channel scope
sc.pp.neighbors(adata, n_neighbors=${sub_neighbors}, n_pcs=${sub_pcs_umap})

# UMAP & t-SNE
print("Computing UMAP...")
sc.tl.umap(adata)
print("Computing t-SNE...")
sc.tl.tsne(adata, use_rep='X_pca')

print("Performing Leiden clustering...")
# Swapped from global params to local input channel scope
sc.tl.leiden(adata, flavor="igraph", resolution=${sub_resolution}, n_iterations=2)

# Save cluster sizes dynamically per sample
cluster_df = adata.obs['leiden'].value_counts().reset_index()
cluster_df.columns = ['cluster', 'n_cells']
cluster_df = cluster_df.sort_values('cluster')
cluster_df.to_csv(f'{file_prefix}_cluster_sizes.csv', index=False)

print(f"Identified {adata.obs['leiden'].nunique()} clusters for test1")
print(cluster_df)

# ==========================================
# 4. CELL TYPE ANNOTATION (RE-RUN FOR SUB-CLUSTERS)
# ==========================================
print("Annotating cell types...")
# Swapped from global params to local input channel scope
gs_list = gene_sets_prepare('${marker_file}', cell_type='${tissue_name}')

sc.pp.scale(adata, max_value=10)
scaled_data = pd.DataFrame(adata.X)
scaled_data.columns = adata.var_names
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

# Save output matches nextflow setup expectations
adata.write(f'{file_prefix}_adata_clustered.h5ad')
sctype_scores.to_csv(f'{file_prefix}_sctype_scores.csv', index=False)
print(f"Cell type annotation complete for test1")
