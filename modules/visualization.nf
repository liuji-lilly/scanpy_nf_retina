process VISUALIZATION {
    tag "Generate visualizations: ${meta.id}"
    publishDir "${params.outdir}/07_visualization", mode: 'copy'
    
    input:
    tuple val(meta), path(adata)
    val(annotated)
    path marker_genes_file
    
    output:
    tuple val(meta), path("*.pdf"), emit: plots
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    import os
    
    # 1. Dynamically append your Nextflow workflow bin path
    sys.path.append(os.path.join("${workflow.projectDir}", "bin"))
    
    # 2. Import your newly written modular python function
    from run_visualization import run_visualization
    
    # 3. Run the visualization securely with typecasting and variable protection
    run_visualization(
        adata_path         = "${adata}",
        sample_id          = "${meta.id}",
        marker_genes_file  = "${marker_genes_file}",
        annotated_str      = "${annotated}"
    )
    """
}