#!/bin/bash

# Example script to run the scanpy pipeline
nextflow run main.nf --input ./sampledata/test --outdir ./sampledata/results --tissue_name "Mouse Organoid"

# Set input directory (10X format)
INPUT_DIR="/home/workspace/scanpytest"

# Set output directory
OUTPUT_DIR="/home/workspace/scanpytest/nextflow_results"

# Run the pipeline with custom parameters
nextflow run main.nf \
    --input ${INPUT_DIR} \
    --outdir ${OUTPUT_DIR} \
    --min_genes 200 \
    --max_genes 2500 \
    --max_mt_percent 5 \
    --n_top_genes 2000 \
    --leiden_resolution 1.0 \
    --annotate_celltypes true \
    -resume \
    -with-report ${OUTPUT_DIR}/pipeline_report.html \
    -with-timeline ${OUTPUT_DIR}/timeline.html \
    -with-dag ${OUTPUT_DIR}/flowchart.svg


