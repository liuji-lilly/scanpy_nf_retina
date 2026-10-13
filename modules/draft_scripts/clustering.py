#!/usr/bin/env python3
import scanpy as sc
import pandas as pd

sc.settings.verbosity = 3

# Load data using the dynamic metadata map pointer
print(f"Loading dimensionality reduction data for ${meta.id}...")
adata = sc.read_h5ad('${adata_dimred}')
#adata=sc.read_h5ad("/fsx/home/l128405/scanpy_nf_retina/sampledata/results/03_dimred/adata_dimred.h5ad")

# Leiden clustering
print("Performing Leiden clustering...")
#sc.tl.leiden(adata, flavor="igraph", resolution=1,n_iterations = 2)
sc.tl.leiden(adata, flavor="igraph", resolution=${resolution},n_iterations = 2)

# Save cluster sizes dynamically per sample
cluster_df = adata.obs['leiden'].value_counts().reset_index()
cluster_df.columns = ['cluster', 'n_cells']
cluster_df = cluster_df.sort_values('cluster')
cluster_df.to_csv('${meta.id}_cluster_sizes.csv', index=False)

print(f"Identified {adata.obs['leiden'].nunique()} clusters for ${meta.id}")
print(cluster_df)

# Save with unique dynamic naming
#adata.write("/fsx/home/l128405/scanpy_nf_retina/sampledata/results/04_clustering/adata_clustered.h5ad")
adata.write('${meta.id}_adata_clustered.h5ad')
