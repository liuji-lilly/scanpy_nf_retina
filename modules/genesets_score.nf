process GENESETS_SCORE {
   tag "pathway & genesets score calculation: ${meta.id}"
    publishDir "${params.outdir}/10_genesets_score", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_annotated)
    path celltype_markergenes_file
    val organism
    val pathway_names
    
    output:
    
    tuple val(meta), path("*.pdf"),          emit: visualization_plots
    
    script:

    """
    #!/usr/bin/env python3
    import sys
    import os
    
    # 1. Dynamically append your Nextflow workflow bin path
    sys.path.append(os.path.join("${workflow.projectDir}", "bin"))
    
    # 2. Import your newly written modular python function
    from run_genesets_score import run_genesets_score
    
    # 3. Run the visualization securely with typecasting and variable protection
    run_genesets_score(
        adata_path         = "${adata_annotated}",
        sample_id          = "${meta.id}",
        organism           = "${organism}",
        marker_file        = "${celltype_markergenes_file}",
        raw_pathway_names  = "${pathway_names}"
    )

    """
}
