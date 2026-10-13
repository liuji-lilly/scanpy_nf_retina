 #!/usr/bin/env python3
import sys
import os

# Dynamically add the directory containing sctypepy to the Python Path
sys.path.append("./bin")

import scanpy as sc
import pandas as pd
from sctypepy.sctype_py import gene_sets_prepare, sctype_score, process_cluster
sc.settings.verbosity = 3

# Load data
adata = sc.read_h5ad('/fsx/home/l128405/scanpy_nf_retina/sampledata/results/05_markers/adata_with_markers.h5ad')

# Prepare gene sets from custom markers file
gs_list = gene_sets_prepare('/fsx/home/l128405/scanpy_nf_retina/bin/sctypepy/ScTypeDB_Full.xlsx', cell_type='Mouse Organoid')

# Scale the data before scoring
sc.pp.scale(adata,max_value=10)
scaled_data = pd.DataFrame(adata.X)
# change column indexes
scaled_data.columns =adata.var_names
# Change the row indexes
scaled_data.index = adata.obs_names
scRNAseqData=scaled_data.T

# Calculate ScType scores for every individual cell
es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])
unique_clusters = adata.obs['leiden'].unique()

# Aggregate scores at the cluster level
# Apply the function to each unique cluster and combine the results into a DataFrame
cL_results = pd.concat([process_cluster(cluster,adata,es_max,'leiden') for cluster in unique_clusters])

# Group by cluster and select the top row based on scores
sctype_scores = cL_results.sort_values('scores', ascending=False).drop_duplicates(subset=['cluster']).sort_values('type').reset_index(drop=True)

# Set low-confidence clusters to "Unknown"
sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

# Iterate over unique clusters
adata.obs['sctype_classification'] = ""
for cluster in sctype_scores['cluster'].unique():
    # Filter sctype_scores for the current cluster
    cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
    # Get the type for the current cluster
    cl_type_value = cl_type['type'].iloc[0]
    # Update 'sctype_classification' in adata.obs for cells belonging to the current cluster
    adata.obs.loc[adata.obs['leiden'] == cluster, 'sctype_classification'] = cl_type_value

# Save
adata.write('adata_annotated.h5ad')
sctype_scores.to_csv('sctype_scores.csv', index=False)
print("Cell type annotation complete")
