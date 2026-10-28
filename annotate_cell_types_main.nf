#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Standalone Annotate Cell Types Runner
========================================================================================
*/
params.tissue_name                  = null
params.marker_file           = null
params.outdir = null
// FIX: Import the SUBCLUSTERING process so Nextflow recognizes it
include { ANNOTATE_CELL_TYPES } from '/fsx/home/l128405/scanpy_nf_retina/modules/annotate.nf'

workflow {
    // 1. Build the primary Data Channel
    final_adata_ch = Channel.fromPath("${params.sub_input_dir}/*_adata_with_markers.h5ad", checkIfExists: true)
        .map { file -> 
            def sample_name = file.baseName.replace('_adata_with_markers', '')
            def meta = [id: sample_name] 
            return tuple(meta, file)
        }

    // 2. Safely read input target CSV and reference marker files
    marker_file_ch = Channel.value(file(params.marker_file, checkIfExists: true))
    tissue_name_val = Channel.value(params.tissue_name)

    ANNOTATE_CELL_TYPES(
        final_adata_ch,
        marker_file_ch,
        tissue_name_val
    )
    final_adata = ANNOTATE_CELL_TYPES.out.adata

}
