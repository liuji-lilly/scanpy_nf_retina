process DIMENSIONALITY_REDUCTION {
    tag "PCA and UMAP: ${meta.id}"
    publishDir "${params.outdir}/03_dimred", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_normalized)
    val(n_pcs)
    val(n_neighbors)
    val(n_pcs_umap)
    
    output:
    tuple val(meta), path("${meta.id}_adata_dimred.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_pca_variance.pdf"), emit: pca_plot
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    sc.settings.verbosity = 3
    # Load data
    print(f"Loading normalized data for ${meta.id}...")
    adata = sc.read_h5ad('${adata_normalized}')
    adata.X = adata.layers["scran_normalization"].copy()
   
    # PCA
    print("Computing PCA...")
    adata.var["highly_variable"] = adata.var["highly_deviant"].copy()
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable", n_comps=${n_pcs})

    # Plot variance explained
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=${n_pcs}, show=False)
    fig = plt.gcf()
    fig.suptitle(f"PCA Variance Ratio - ${meta.id}", y=1.02, fontsize=12)
    plt.savefig('${meta.id}_pca_variance.pdf',format="pdf", bbox_inches="tight", dpi=300)
    plt.close()

    # Compute neighborhood graph
    print("Computing neighborhood graph...")
    ##n_neighbors = 15
    #n_pcs_umap = 40
    sc.pp.neighbors(adata, n_neighbors=${n_neighbors}, n_pcs=${n_pcs_umap})


    # UMAP
    print("Computing UMAP...")
    sc.tl.umap(adata)
    print("Computing t-SNE...")
    sc.tl.tsne(adata, use_rep='X_pca')

    # Save
    #adata.write("/fsx/home/l128405/scanpy_nf_retina/sampledata/results/03_dimred/adata_normalized.h5ad")
    adata.write('${meta.id}_adata_dimred.h5ad')
    print(f"Dimensionality reduction complete for ${meta.id}")

    """
}