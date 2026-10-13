#!/usr/bin/env python3
import scanpy as sc
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import logging
import matplotlib.pyplot as plt
import numpy as np
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import seaborn as sns
from rpy2.robjects import default_converter, numpy2ri, pandas2ri, r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy.sparse import issparse

scry = importr('scry')
methods = importr('methods')
matrix = importr('Matrix')
single_cell_experiment = importr('SingleCellExperiment')
seurat = importr('Seurat')
scater = importr('scater')
scran = importr('scran')
bioc_parallel = importr('BiocParallel')


sc.settings.verbosity = 3
# Set figure parameters for clean, minimal plots
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)


# Load data using the dynamic metadata dictionary pointer
print(f"Loading QC data for ${meta.id}...")
adata = sc.read_h5ad('${adata_qc}')
#adata = sc.read_h5ad("/fsx/home/l128405/scanpy_nf_retina/sampledata/results/01_qc/test1_qc.h5ad")
adata.layers["counts"] = adata.X.copy()
# Normalize

print("Normalizing data by scran...")
#scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
#adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(
    adata_pp, key_added="groups", flavor="igraph", n_iterations=2, directed=False
)
data_mat = adata.layers["counts"].T
# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
    # Convert sparse matrix to dense numpy array
    data_mat = data_mat.toarray()

with localconverter(ro.default_converter + numpy2ri.converter):
    ro.globalenv["data_mat"] = data_mat

with localconverter(ro.default_converter + pandas2ri.converter):
    ro.globalenv["input_groups"] = adata_pp.obs["groups"]

del adata_pp

ro.r('''

size_factors = sizeFactors(
    computeSumFactors(
        SingleCellExperiment(
            list(counts=data_mat)), 
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
    )
)

''')
adata.obs["size_factors"]= np.array(ro.globalenv["size_factors"])
scran = adata.X / adata.obs["size_factors"].values[:, None]
scran_logged = np.log1p(scran)
adata.layers["scran_normalization"] = csr_matrix(scran_logged)

#print("Normalizing data by pearson residuals...")
analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

# Identify highly variable genes
print("Identifying highly variable genes...")
# # Initialize R packages directly into Python variables

X_sparse = adata.X.T.tocoo()

with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    ro.globalenv["obs"] = adata.obs
    ro.globalenv["var"] = adata.var

i, j = X_sparse.row, X_sparse.col
x = X_sparse.data

ro.globalenv["i"] = ro.IntVector((i + 1).tolist())  # R is 1-indexed
ro.globalenv["j"] = ro.IntVector((j + 1).tolist())
ro.globalenv["x"] = ro.FloatVector(x.tolist())

r("X <- sparseMatrix(i = i, j = j, x = x, dims = c({}, {}))".format(*X_sparse.shape))

ro.r('''

sce <- SingleCellExperiment(
  assays = list(X = X),
  colData = obs,
  rowData = var
)

sce <- devianceFeatureSelection(sce, assay = "X")

''')

with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
    binomial_deviance = ro.r("rowData(sce)$binomial_deviance")

idx = binomial_deviance.argsort()[-${n_top_genes}:]
#idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


sc.pp.highly_variable_genes(adata, layer="scran_normalization")

plt.figure(figsize=(6, 5))
ax = sns.scatterplot(
    data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=15
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
# Save the figure to a PDF file
plt.savefig('${meta.id}_hvg_plot.pdf', format="pdf", bbox_inches="tight", dpi=300)
plt.close()


## Regress out and scale
#print("Regressing out effects...")
#sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#print("Scaling data...")
#sc.pp.scale(adata, max_value=10)

# Save with unique dynamic naming
#adata.write("/fsx/home/l128405/scanpy_nf_retina/sampledata/results/02_normalize/adata_normalized.h5ad")
import anndata
# 1. Turn on the setting to allow writing these string types
anndata.settings.allow_write_nullable_strings = True
adata.write('${meta.id}_adata_normalized.h5ad')
print(f"Normalization complete for ${meta.id}: {adata.n_obs} cells, {adata.n_vars} genes")
