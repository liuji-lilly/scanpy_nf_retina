process NORMALIZE_AND_HVG {
    tag "Normalization and HVG: ${meta.id}"
    publishDir "${params.outdir}/02_normalize", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_qc)
    val(n_top_genes)
    
    output:
    tuple val(meta), path("${meta.id}_adata_normalized.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_hvg_plot.pdf"),          emit: hvg_plot
    
    script:
    """
    #!/usr/bin/env python3
    import matplotlib
    matplotlib.use('Agg')  # FIX: Moved to the absolute top before pyplot initialization
    import matplotlib.pyplot as plt
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import logging
    import seaborn as sns
    from scipy.sparse import csr_matrix, issparse
    
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects import default_converter, numpy2ri, pandas2ri, r
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    # Initialize R packages directly into runtime environments
    scry = importr('scry')
    methods = importr('methods')
    matrix = importr('Matrix')
    single_cell_experiment = importr('SingleCellExperiment')
    seurat = importr('Seurat')
    scater = importr('scater')
    scran = importr('scran')
    bioc_parallel = importr('BiocParallel')

    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

    print("Loading QC data for ${meta.id}...")
    adata = sc.read_h5ad('${adata_qc}')
    adata.layers["counts"] = adata.X.copy()
    
    # Pre-cluster cells briefly to improve Scran pool factor stability
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups", flavor="igraph", n_iterations=2, directed=False)
    
    data_mat = adata.layers["counts"].T
    # Optimize matrix data layout tracking styles to resolve memory inflation anomalies
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
        data_mat = data_mat.toarray()

    with localconverter(ro.default_converter + numpy2ri.converter):
        ro.globalenv["data_mat"] = data_mat

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv["input_groups"] = adata_pp.obs["groups"]

    del adata_pp

    ro.r('''
    library(BiocParallel)
    size_factors = sizeFactors(
        computeSumFactors(
            SingleCellExperiment(list(counts=data_mat)), 
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
        )
    )
    ''')
    
    adata.obs["size_factors"] = np.array(ro.globalenv["size_factors"])
    scran_mat = adata.X / adata.obs["size_factors"].values[:, None]
    scran_logged = np.log1p(scran_mat)
    adata.layers["scran_normalization"] = csr_matrix(scran_logged)

    # Compute Pearson Residuals secondary matrix layers
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

    # ====================================================================
    # BINOMIAL DEVIANCE FEATURE SELECTION (SCRY)
    # ====================================================================
    print("Identifying highly deviant features via scry...")
    X_sparse = adata.X.T.tocoo()

    with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
        ro.globalenv["obs"] = adata.obs
        ro.globalenv["var"] = adata.var

    i, j = X_sparse.row, X_sparse.col
    x = X_sparse.data

    ro.globalenv["i"] = ro.IntVector((i + 1).tolist()) 
    ro.globalenv["j"] = ro.IntVector((j + 1).tolist())
    ro.globalenv["x"] = ro.FloatVector(x.tolist())

    # FIX: Swapped .format() curly brackets to old-school python modulo syntax to prevent Nextflow from throwing template compilation crashes
    r("X <- sparseMatrix(i = i, j = j, x = x, dims = c(%d, %d))" % X_sparse.shape)

    ro.r('''
    sce <- SingleCellExperiment(assays = list(X = X),colData = obs,rowData = var)
    sce <- devianceFeatureSelection(sce, assay = "X")
    ''')

    with localconverter(default_converter + pandas2ri.converter + numpy2ri.converter):
        binomial_deviance = ro.r("rowData(sce)\$binomial_deviance")

    idx = binomial_deviance.argsort()[-${n_top_genes}:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True

    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance

    # Run native highly variable gene checking on your new layer
    sc.pp.highly_variable_genes(adata, layer="scran_normalization")
    
    ##sc.pp.highly_variable_genes(adata, n_top_genes=${n_top_genes}, flavor="seurat_v3", layer="counts")
    plt.figure(figsize=(6, 5))
    ax = sns.scatterplot(data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=15)
    ax.set_xlim(None, 1.5)
    ax.set_ylim(None, 3)
    plt.savefig('${meta.id}_hvg_plot.pdf', format="pdf", bbox_inches="tight", dpi=300)
    plt.close()

    import anndata
    anndata.settings.allow_write_nullable_strings = True
    adata.write('${meta.id}_adata_normalized.h5ad')
    print(f"Normalization complete for ${meta.id}: {adata.n_obs} cells")
    """
}