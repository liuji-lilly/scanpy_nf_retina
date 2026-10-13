process GENERATE_ALL_REPORTS {
    tag "Generate master summary report"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path("adatas/*")
    path("qc_metrics/*")
    path("markers/*")
    
    output:
    path("combined_analysis_summary.html"), emit: report
    path("combined_analysis_summary.md"),   emit: markdown
    
    script:
    """
    # Your python script will need to be updated to look inside these directories
    report_generator.py --adata_dir adatas --qc_dir qc_metrics --marker_dir markers
    """
}