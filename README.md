# Scanpy Single-Cell RNA-seq Nextflow Pipeline

A comprehensive, production-ready Nextflow pipeline for single-cell RNA-seq analysis using Scanpy.

**Repository**: https://github.com/liuji-lilly/scanpy_nf_retina

## Installation

### Quick Start
```bash
# Clone the repository
git clone https://github.com/liuji-lilly/scanpy_nf_retina.git
cd scanpy_nf_retina

# Run with conda environment (recommended)
nextflow run main.nf -profile conda --input /path/to/10x_data --outdir results
```

### Prerequisites
- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- One of:
  - Conda/Mamba 
  - Docker
  - Singularity
  - Python 3.10+ with scanpy installed

## Features

- Fully automated workflow from raw 10X data to annotated cell types
- Quality control with customizable thresholds
- Normalization and highly variable gene selection
- PCA and UMAP dimensionality reduction
- Leiden clustering
- Marker gene identification of each cluster
- Cell type annotation
- Comprehensive visualizations 
- HTML and Markdown reports

## Environment Management

The pipeline supports **4 execution modes** to match your needs:

1. **Standard** (default) - Uses your current Python environment
   ```bash
   nextflow run main.nf --input data --outdir results
   ```

2. **Conda** (recommended) - Nextflow manages conda environment automatically
   ```bash
   nextflow run main.nf -profile conda --input data --outdir results
   ```

3. **Docker** - Maximum reproducibility with containers
   ```bash
   nextflow run main.nf -profile docker --input data --outdir results
   ```

4. **Singularity** - For HPC clusters
   ```bash
   nextflow run main.nf -profile singularity --input data --outdir results
   ```

📖 **See [ENVIRONMENT.md](ENVIRONMENT.md) for complete setup instructions** including:
- Detailed installation steps for each mode
- Conda environment.yml
- Dockerfile and build scripts
- HPC cluster configurations
- Troubleshooting guide

## Requirements

### Software
- Nextflow >= 23.04.0
- Python >= 3.9
- Scanpy >= 1.9.0

### Python packages
```bash
pip install scanpy matplotlib pandas openpyxl
pip install igraph leidenalg  # For clustering
```

## Quick Start

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
./nextflow -version
```

### 2. Run the pipeline

```bash
# Basic usage with default parameters
nextflow run main.nf --input /path/to/10x/data --outdir results

# With custom parameters
nextflow run main.nf --input ./sampledata/samplesheet.csv --outdir ./sampledata/results --tissue_name "Mouse Organoid"  
```

### 3. Resume a failed run

```bash
nextflow run main.nf --input /path/to/10x/data --outdir results -resume
```

## Input Data

The pipeline expects 10X Genomics format data:
```
input_directory/
├── matrix.mtx      # Expression matrix
├── genes.tsv       # Gene information
└── barcodes.tsv    # Cell barcodes
```

## Parameters

### Input/Output
- `--input`: sample sheet
Format:
```
sample_id,file_path,batch,treatment,chemistry
retina_P14,/fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment/cellranger/count/retina_P14/outs/filtered_feature_bc_matrix,Li_2024,No enrichment,Next GEM
retina_P90,/fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment/cellranger/count/retina_P90/outs/filtered_feature_bc_matrix,Li_2024,No enrichment,Next GEM
retina_P120,/fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment2/cellranger/count/retina_P120/outs/filtered_feature_bc_matrix,Li_2024,No enrichment,Next GEM
retina_P40,/fsx/home/l128405/10xgenomics/retina_dataset/cellranger_alignment2/cellranger/count/retina_P40/outs/filtered_feature_bc_matrix,Li_2024,No enrichment,Next GEM
```
- `--outdir`: Output directory (default: `./results`)

### QC Parameters
- `--min_genes`: Minimum genes per cell (default: 200)
- `--min_cells`: Minimum cells per gene (default: 3)
- `--max_genes`: Maximum genes per cell (default: 2500)
- `--max_mt_percent`: Maximum mitochondrial percentage (default: 5)

### Analysis Parameters
- `--n_top_genes`: Number of highly variable genes (default: 2000)
- `--n_pcs`: Number of principal components (default: 50)
- `--n_neighbors`: Number of neighbors for graph (default: 10)
- `--n_pcs_umap`: PCs used for UMAP (default: 40)
- `--leiden_resolution`: Leiden clustering resolution (default: 1.0)
- `--annotate_celltypes`: Enable cell type annotation (default: true)

### Resource Parameters
- `--max_cpus`: Maximum CPUs (default: 4)
- `--max_memory`: Maximum memory (default: '16.GB')
- `--max_time`: Maximum time (default: '4.h')

## Output Structure

```
results/
├── 01_qc
│   ├── retina_P120_qc.h5ad
│   ├── retina_P120_qc_cutoffs.csv
│   ├── retina_P120_qc_metrics.csv
│   ├── retina_P120_qc_plots.pdf
│   ├── retina_P120_raw.h5ad
├── 02_normalize
│   ├── retina_P120_adata_normalized.h5ad
│   ├── retina_P120_hvg_plot.pdf
├── 03_dimred
│   ├── retina_P120_adata_dimred.h5ad
│   ├── retina_P120_pca_variance.pdf
├── 04_clustering
│   ├── retina_P120_adata_clustered.h5ad
│   ├── retina_P120_cluster_sizes.csv
├── 05_markers
│   ├── retina_P120_adata_with_markers.h5ad
│   ├── retina_P120_marker_genes.csv
│   ├── retina_P120_marker_genes_dotplot.png
│   ├── retina_P120_marker_genes_heatmap.png
├── 06_annotation
│   ├── pipeline_info
│   │   ├── execution_trace.txt
│   │   └── pipeline_dag.dot
│   ├── retina_P120_adata_annotated.h5ad
│   ├── retina_P120_cell_type_summary.csv
│   ├── retina_P120_sctype_scores.csv
├── 07_visualization
│   ├── retina_P120_adata_all_markers_individual_violins.pdf
│   ├── retina_P120_adata_global_celltype_composition.pdf
│   ├── retina_P120_adata_marker_custom_pca_expression.pdf
│   ├── retina_P120_adata_marker_custom_tsne_expression.pdf
│   ├── retina_P120_adata_marker_custom_umap_expression.pdf
│   ├── retina_P120_adata_marker_selectedmarkers_dotplot.pdf
│   ├── retina_P120_adata_marker_selectedmarkers_dotplot_scaled.pdf
│   ├── retina_P120_adata_marker_selectedmarkers_stacked_violin.pdf
│   ├── retina_P120_adata_pca_clusters.pdf
│   ├── retina_P120_adata_pca_comprehensive_qc.pdf
│   ├── retina_P120_adata_tsne_clusters.pdf
│   ├── retina_P120_adata_tsne_comprehensive_qc.pdf
│   ├── retina_P120_adata_umap_cell_types.pdf
│   ├── retina_P120_adata_umap_comprehensive_qc.pdf
│   ├── retina_P120_adata_violin_comprehensive_qc.pdf
├── 10_genesets_score
│   ├── retina_P120_adata_celltype_aucell_individual_stripplots.pdf
│   ├── retina_P120_adata_celltype_aucell_individual_violins.pdf
│   ├── retina_P120_adata_progeny_all_pathways_dotplot.pdf
│   ├── retina_P120_adata_progeny_all_pathways_pca.pdf
│   ├── retina_P120_adata_progeny_all_pathways_stacked_violin.pdf
│   ├── retina_P120_adata_selected_pathway_dotplot.pdf
│   ├── retina_P120_adata_selected_pathway_individual_violins.pdf
│   ├── retina_P120_adata_selected_pathway_pca.pdf
│   ├── retina_P120_adata_selected_pathway_stacked_violin.pdf
├── combined_analysis_summary.html
├── combined_analysis_summary.md
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   └── pipeline_dag.dot
```

## Execution Profiles

### Local (default)
```bash
nextflow run main.nf --input data --outdir results
```

### Docker
```bash
nextflow run main.nf --input data --outdir results -profile docker
```

### Singularity
```bash
nextflow run main.nf --input data --outdir results -profile singularity
```

### SLURM
```bash
nextflow run main.nf --input data --outdir results -profile slurm
```

Base your annotations on the marker genes identified in `05_markers/marker_genes.csv`.

## Advanced Usage

### Running specific steps only

You can modify the `main.nf` to skip certain steps by commenting them out.

### Adjusting cluster resolution

Higher resolution = more clusters:
```bash
nextflow run main.nf --input data --leiden_resolution 1.5
```

Lower resolution = fewer clusters:
```bash
nextflow run main.nf --input data --leiden_resolution 0.5
```

### Parallel execution on HPC

Create a custom profile in `nextflow.config`:

```groovy
profiles {
    hpc {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.cpus = 8
        process.memory = '32.GB'
        process.time = '8.h'
    }
}
```

Run with: `nextflow run main.nf -profile hpc`

## Troubleshooting

### Out of memory errors
Increase memory: `--max_memory '32.GB'`

### Clustering fails
Install igraph: `pip install igraph leidenalg`

### Pipeline hangs
Check logs in `work/` directory

### Resume from checkpoint
Always use `-resume` flag to continue from where it stopped

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology.
- **Scanpy**: Wolf, F.A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology.

## Support

For issues and questions:
- Check the [Scanpy documentation](https://scanpy.readthedocs.io/)
- Check the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

