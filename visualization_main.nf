#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Standalone Visualization Runner
========================================================================================
*/

// FIX: Import the VISUALIZATION process so Nextflow recognizes it
include { VISUALIZATION } from './modules/visualization.nf'
params.outdir = null
params.sub_input_dir      = null
params.custom_marker_file = null
params.annotate_celltypes = true


workflow {
    // FIX 1: Provide absolute system path targeting your actual single-cell data file
    // Adjust the path inside the quotes if your file isn't directly inside 'results/06_annotation/'
    // 1. Build the primary Data Channel
    final_adata_ch = Channel.fromPath("${params.sub_input_dir}/*_annotated.h5ad", checkIfExists: true)
        .map { file -> 
            def sample_name = file.baseName.replace('_annotated', '')
            def meta = [id: sample_name] 
            return tuple(meta, file)
        }

    custom_markers_ch = Channel.value(file(params.custom_marker_file, checkIfExists: true))

    VISUALIZATION(
        final_adata_ch,
        params.annotate_celltypes,
        custom_markers_ch
    )
}