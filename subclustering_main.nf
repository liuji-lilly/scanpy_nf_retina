#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Standalone Subclustering Runner
========================================================================================
*/

// FIX: Import the SUBCLUSTERING process so Nextflow recognizes it
include { SUBCLUSTERING } from './modules/subclustering.nf'

workflow {
    // 1. Build the primary Data Channel
    final_adata_ch = Channel.fromPath("${params.sub_input_dir}/*_annotated.h5ad", checkIfExists: true)
        .map { file -> 
            def sample_name = file.baseName.replace('_annotated', '')
            def meta = [id: sample_name] 
            return tuple(meta, file)
        }

    // 2. Safely read input target CSV and reference marker files
    ch_custom_markers     = params.custom_marker_file ? file(params.custom_marker_file, checkIfExists: true) : []
    ch_marker_file      = params.marker_file ? file(params.marker_file, checkIfExists: true) : []
    ch_interest_csv         = params.interest_csv ? file(params.interest_csv, checkIfExists: true) : []
    // 3. Clean signature execution call (No positional mapping risks)
    SUBCLUSTERING(
        final_adata_ch,
        ch_interest_csv,
        ch_custom_markers,
        ch_marker_file,
        params.annotate_celltypes,
        params.tissue_name
    )
}
