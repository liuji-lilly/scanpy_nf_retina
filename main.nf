#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Single-cell RNA-seq Analysis Pipeline using Scanpy
========================================================================================
    Author: Xiaobin Zheng
    Description: Automated pipeline for scRNA-seq QC, normalization, clustering, 
                 and cell type identification
========================================================================================
*/

// --- PIPELINE PARAMETER DEFAULTS ---
params.input                        = null
params.outdir                       = "./results"
params.n_top_genes                  = 2000
params.n_pcs                        = 30
params.n_neighbors                  = 15
params.n_pcs_umap                   = 30
params.leiden_resolution            = 1.0
params.annotate_celltypes           = true
params.marker_file                  = null
params.tissue_name                  = null
params.custom_marker_file           = null

// Labeled bindings matching your advanced MAD filtering overrides
params.min_cells                    = 3
params.min_genes                    = 200
params.max_genes                    = 3000
params.min_transcripts              = 500
params.max_transcripts              = 10000
params.max_pct_counts_in_top_20_genes = 80.0
params.max_mt_percent               = 5.0




// Print pipeline header
def printHeader() {
    log.info """\\
    ===================================================
     S C A N P Y   P I P E L I N E
    ===================================================
     Input directory : ${params.input}
     Output directory: ${params.outdir}
     Min cells/gene  : ${params.min_cells}
     Min genes/cell  : ${params.min_genes}
     Max genes/cell  : ${params.max_genes}
     Min transcripts : ${params.min_transcripts}
     Max transcripts : ${params.max_transcripts}
     Max Top 20%     : ${params.max_pct_counts_in_top_20_genes}
     Max MT%         : ${params.max_mt_percent}
     Leiden res.     : ${params.leiden_resolution}
     Annotate cells  : ${params.annotate_celltypes}
     Marker file     : ${params.marker_file}
     Tissue name     : ${params.tissue_name}
     Custom plot file: ${params.custom_marker_file}
    ===================================================
     """.stripIndent()
}

// Validate parameters
def validateParams() {
    if (!params.input) {
        error "Please provide an input directory samplesheet via --input"
    }
    if (!file(params.input).exists()) {
        error "Input samplesheet does not exist: ${params.input}"
    }
    if (params.annotate_celltypes) {
        if (!params.marker_file) {
            error "Cell type annotation is enabled. Please provide a marker database file via --marker_file"
        }
        if (!file(params.marker_file).exists()) {
            error "Provided marker file does not exist: ${params.marker_file}"
        }
        if (!params.tissue_name) {
            error "Cell type annotation is enabled. Please specify a tissue option via --tissue_name (e.g., 'Liver', 'Eye', 'Brain')"
        }
    }
    if (!params.custom_marker_file) {
        error "Please provide a custom marker genes CSV file via --custom_marker_file for visualization"
    }
    if (!file(params.custom_marker_file).exists()) {
        error "Provided custom marker file does not exist: ${params.custom_marker_file}"
    }
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { LOAD_AND_QC } from './modules/load_qc.nf'
include { NORMALIZE_AND_HVG } from './modules/normalize.nf'
include { DIMENSIONALITY_REDUCTION } from './modules/dimensionality_reduction.nf'
include { CLUSTERING } from './modules/clustering.nf'
include { MARKER_GENES } from './modules/marker_genes.nf'
include { ANNOTATE_CELL_TYPES } from './modules/annotate.nf'
include { VISUALIZATION } from './modules/visualization.nf'
include { GENERATE_ALL_REPORTS } from './modules/report.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/
workflow {
    printHeader()
    validateParams()
    
    // 1. Create a parallelized channel parsing the CSV samplesheet.
    samples_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            def meta = [:]
            meta.id        = row.sample_id
            meta.batch     = row.batch
            meta.treatment = row.treatment
            meta.chemistry = row.chemistry
            def data_dir = file(row.file_path, checkIfExists: true)
            return tuple(meta, data_dir)
        }
        
    // Step 1: Load data and perform MAD filtering + Outlier evaluation
    // FIXED: Passed exactly 8 ordered argument layers mapping to your modular load_qc update
    LOAD_AND_QC(
        samples_ch
    )

    // Step 2: Normalize by Scran pools & pearson metrics + Identify highly variable features via scry
    // FIXED: Targets '.out.adata_qc' channel emmission label matching your new module name tag
    NORMALIZE_AND_HVG(
        LOAD_AND_QC.out.adata_qc,
        params.n_top_genes
    )
    
    // Step 3: Dimensionality reduction (PCA and UMAP)
    DIMENSIONALITY_REDUCTION(
        NORMALIZE_AND_HVG.out.adata,
        params.n_pcs,
        params.n_neighbors,
        params.n_pcs_umap
    )
    
    // Step 4: Clustering
    CLUSTERING(
        DIMENSIONALITY_REDUCTION.out.adata,
        params.leiden_resolution
    )
    
    // Step 5: Identify marker genes
    MARKER_GENES(
        CLUSTERING.out.adata
    )
    
    // Step 6: Annotate cell types (if enabled)
    if (params.annotate_celltypes) {
        marker_file_ch = Channel.value(file(params.marker_file, checkIfExists: true))
        tissue_name_val = Channel.value(params.tissue_name)
        species_val     = Channel.value(params.species)

        ANNOTATE_CELL_TYPES(
            MARKER_GENES.out.adata,
            marker_file_ch,
            tissue_name_val,
            species_val
        )
        final_adata = ANNOTATE_CELL_TYPES.out.adata
    } else {
        final_adata = MARKER_GENES.out.adata
    }
    
    custom_markers_ch = Channel.value(file(params.custom_marker_file, checkIfExists: true))

    // Step 7: Generate visualizations
    VISUALIZATION(
        final_adata,
        params.annotate_celltypes,
        custom_markers_ch
    )

    // Step 8: Generate summary report
    // FIXED: Stripped meta keys and collected separate outputs into clean group arrays
    ch_all_adatas   = final_adata.map { meta, file -> file }.collect()
    ch_all_qc       = LOAD_AND_QC.out.qc_metrics
                        .map { meta, file -> file }
                        .flatten()
                        .filter { it.name.endsWith('_qc_metrics.csv') }
                        .collect()
    ch_all_markers  = MARKER_GENES.out.markers.map { meta, file -> file }.collect()

    GENERATE_ALL_REPORTS(
        ch_all_adatas,
        ch_all_qc,
        ch_all_markers
    )
}