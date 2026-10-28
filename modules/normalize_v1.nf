process NORMALIZE_AND_HVG {
    tag "Normalization and HVG: ${meta.id}"
    publishDir "${params.outdir}/02_normalize", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_qc)
    val(n_top_genes)
    
    output:
    tuple val(meta), path("${meta.id}_adata_normalized.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_hvg_plot.png"), emit: hvg_plot
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    sc.settings.verbosity = 3
    
    # Load data using the dynamic metadata dictionary pointer
    print(f"Loading QC data for ${meta.id}...")
    adata = sc.read_h5ad('${adata_qc}')
    
    # Normalize
    print("Normalizing data...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    print("Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=${n_top_genes})
    
    # Plot HVGs
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.suptitle(f"Highly Variable Genes - ${meta.id}", y=1.02, fontsize=12)
    plt.savefig('${meta.id}_hvg_plot.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    ## Regress out and scale
    #print("Regressing out effects...")
    #sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    #print("Scaling data...")
    #sc.pp.scale(adata, max_value=10)
    
    # Save with unique dynamic naming
    adata.write('${meta.id}_adata_normalized.h5ad')
    print(f"Normalization complete for ${meta.id}: {adata.n_obs} cells, {adata.n_vars} genes")
    """
}
