process DIFFERENTIAL_GENE_EXPRESSION {
    tag "Subsetting & Reclustering: ${meta.id}"
    publishDir "${params.outdir}/08_subclustering", mode: 'copy'
    
    input:
    tuple val(meta), path(adata_annotated)
    path interest_csv
    path custom_marker_file
    path marker_file
    val(annotated)
    val(tissue_name)
    
    output:
    tuple val(meta), path("${meta.id}.*_adata_clustered.h5ad"),  emit: adata
    tuple val(meta), path("${meta.id}.*_cluster_sizes.csv"),     emit: cluster_sizes
    tuple val(meta), path("${meta.id}.*_hvg_plot.pdf"),          emit: hvg_plot
    tuple val(meta), path("${meta.id}.*_pca_variance.pdf"),      emit: pca_plot
    tuple val(meta), path("${meta.id}.*_sctype_scores.csv"),     emit: sctype_scores
    tuple val(meta), path("${meta.id}.*.pdf"),          emit: visualization_plots
    script:
    """

    """
}
