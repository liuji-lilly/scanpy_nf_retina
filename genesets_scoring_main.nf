#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Standalone Subclustering Runner
========================================================================================
*/

// FIX: Import the GENESETS_SCORE process so Nextflow recognizes it
include { GENESETS_SCORE } from '/fsx/home/l128405/scanpy_nf_retina/modules/genesets_score.nf'

workflow {
    // 1. Build the primary Data Channel
    final_adata_ch = Channel.fromPath("${params.sub_input_dir}/*_annotated.h5ad", checkIfExists: true)
        .map { file -> 
            def sample_name = file.baseName.replace('_annotated', '')
            def meta = [id: sample_name] 
            return tuple(meta, file)
        }

    // 2. Safely read input target CSV and reference marker files
    ch_marker_file     = params.celltype_markergenes_file ? file(params.celltype_markergenes_file, checkIfExists: true) : []
    ch_interest_csv         = params.interest_csv ? file(params.interest_csv, checkIfExists: true) : []
    // 3. Clean signature execution call (No positional mapping risks)
    GENESETS_SCORE(
        final_adata_ch,
        ch_marker_file,
        params.organism,
        params.pathway_names
    )
}
