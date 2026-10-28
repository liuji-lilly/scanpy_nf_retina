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
    --n_neighbors 10 \
    --n_top_genes 4000 \
    --leiden_resolution 2.0 \
    -resume \
    -with-report ${OUTPUT_DIR}/pipeline_report.html \
    -with-timeline ${OUTPUT_DIR}/timeline.html \
    -with-dag ${OUTPUT_DIR}/flowchart.svg



# Set up paths
CELLRANGER_REF="/fsx/home/l128405/10xgenomics/refdata-gex-GRCm39-2024-A"
samplesheet="/fsx/home/l128405/10xgenomics/retina_dataset/samplesheet.csv"
outputdir="/fsx/home/l128405/10xgenomics/retina_dataset/results"
nextflow_config="/fsx/home/l128405/10xgenomics/retina_dataset/nextflow.config"

# Execute Nextflow orchestration
# nextflow run nf-core/scrnaseq \
#   -profile singularity \
#   -c $nextflow_config \
#   --input $samplesheet \
#   --cellranger_index $CELLRANGER_REF \
#   --protocol auto \
#   --aligner cellranger \
#   --outdir $outputdir \
#   -resume

CELLRANGER_REF="/fsx/home/l128405/10xgenomics/refdata-gex-GRCm39-2024-A"
nextflow_config="/fsx/home/l128405/10xgenomics/retina_dataset/nextflow_gpu.config"
outputdir="/fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment"
samplesheet="/fsx/home/l128405/10xgenomics/retina_dataset/samplesheet.csv"

mkdir -p /fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment
cd /fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment

  # Execute Nextflow orchestration
nextflow run nf-core/scrnaseq \
  -profile singularity \
  -c $nextflow_config \
  --input $samplesheet \
  --cellranger_index $CELLRANGER_REF \
  --protocol auto \
  --aligner cellranger \
  --outdir $outputdir \
  -resume