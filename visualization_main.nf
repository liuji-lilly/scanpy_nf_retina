#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Standalone Visualization Runner
========================================================================================
*/

// FIX: Import the VISUALIZATION process so Nextflow recognizes it
include { VISUALIZATION } from './modules/visualization.nf'


workflow {
    // FIX 1: Provide absolute system path targeting your actual single-cell data file
    // Adjust the path inside the quotes if your file isn't directly inside 'results/06_annotation/'
    final_adata_ch = Channel.fromPath("${projectDir}/sampledata/results/06_annotation/*.h5ad", checkIfExists: true)
        .map { file -> 
            // This safely extracts 'adata_annotated' or the actual sample name
            def sample_name = file.baseName
            def meta = [id: sample_name] 
            return tuple(meta, file)
        }

    custom_markers_ch = Channel.fromPath(params.custom_marker_file, checkIfExists: true)

    VISUALIZATION(
        final_adata_ch,
        params.annotate_celltypes,
        custom_markers_ch
    )
}