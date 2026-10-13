process CLUSTERING {
    tag "Leiden clustering: ${meta.id}"
    publishDir "${params.outdir}/04_clustering", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_dimred)
    val(resolution)
    
    output:
    tuple val(meta), path("${meta.id}_adata_clustered.h5ad"), emit: adata
    tuple val(meta), path("${meta.id}_cluster_sizes.csv"), emit: cluster_sizes
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    
    sc.settings.verbosity = 3
    
    # Load data using the dynamic metadata map pointer
    print(f"Loading dimensionality reduction data for ${meta.id}...")
    adata = sc.read_h5ad('${adata_dimred}')
    
    # Leiden clustering
    print("Performing Leiden clustering...")
    sc.tl.leiden(adata, resolution=${resolution})
    
    # Save cluster sizes dynamically per sample
    cluster_df = adata.obs['leiden'].value_counts().reset_index()
    cluster_df.columns = ['cluster', 'n_cells']
    cluster_df = cluster_df.sort_values('cluster')
    cluster_df.to_csv('${meta.id}_cluster_sizes.csv', index=False)
    
    print(f"Identified {adata.obs['leiden'].nunique()} clusters for ${meta.id}")
    print(cluster_df)
    
    # Save with unique dynamic naming
    adata.write('${meta.id}_adata_clustered.h5ad')
    """
}