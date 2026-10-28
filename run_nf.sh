conda activate scanpy-pipeline 

nextflow run main.nf --help

nextflow run main.nf \
--input ./sampledata/samplesheet.csv \
--outdir ./sampledata/results --tissue_name "Mouse Organoid" \
--custom_marker_file "/fsx/home/l128405/scanpy_nf_retina/sampledata/marker_genes.csv"

nextflow run main.nf \
--input ./sampledata/samplesheet.csv \
--outdir ./sampledata/results --tissue_name "Mouse Organoid" \
-resume 

cd /fsx/home/l128405/10xgenomics/multiome/S41A/
samplesheet=/fsx/home/l128405/10xgenomics/multiome/S41A/processing_samplesheet.csv
outputdir=/fsx/home/l128405/10xgenomics/multiome/S41A/analysis
markergene_file=/fsx/home/l128405/scanpy_nf_retina/bin/marker_genes.csv
customconfig=/fsx/home/l128405/scanpy_nf_retina/nextflow.config

nextflow run /fsx/home/l128405/scanpy_nf_retina/main.nf \
--input $samplesheet \
-c $customconfig \
--outdir $outputdir --tissue_name "P20 sc-Cochlea" \
--custom_marker_file $markergene_file \
--n_top_genes 4000 \
--n_neighbors 10 \
--leiden_resolution 2.0 \
-resume

nextflow run visualization_main.nf \
--custom_marker_file /fsx/home/l128405/scanpy_nf_retina/sampledata/marker_genes.csv   \
--annotate_celltypes false   --outdir "./results"

nextflow run subclustering_main.nf \
  -c nextflow.config \
  --interest_csv "/fsx/home/l128405/scanpy_nf_retina/bin/subclustering.csv" \
  --sub_input_dir "/fsx/home/l128405/scanpy_nf_retina/sampledata/results/06_annotation/" \
  --custom_marker_file "/fsx/home/l128405/scanpy_nf_retina/sampledata/marker_genes.csv" \
  --tissue_name "Mouse Organoid"

cd /fsx/home/l128405/10xgenomics/retina_dataset/Li_2024
nextflow run /fsx/home/l128405/scanpy_nf_retina/annotate_cell_types_main.nf \
  -c /fsx/home/l128405/scanpy_nf_retina/nextflow.config \
  --sub_input_dir "/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024/05_markers" \
  --outdir "/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024" \
  --marker_file "/fsx/home/l128405/scanpy_nf_retina/bin/sctypepy/ScTypeDB_Full.xlsx" \
  --tissue_name "Eye"

custom_marker_file="/fsx/home/l128405/10xgenomics/retina_dataset/retina_markergenes.csv"
nextflow run /fsx/home/l128405/scanpy_nf_retina/visualization_main.nf \
--sub_input_dir "/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024/06_annotation" \
--custom_marker_file $custom_marker_file   \
--annotate_celltypes false \
--outdir "/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024"

nextflow run /fsx/home/l128405/scanpy_nf_retina/genesets_scoring_main.nf \
  -c /fsx/home/l128405/scanpy_nf_retina/nextflow.config \
  --sub_input_dir "/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024/06_annotation" \
  --outdir "/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024" \
  --celltype_markergenes_file "/fsx/home/l128405/10xgenomics/retina_dataset/retina_markergenes.csv" \
  --organism "mouse" --pathway_names "Wnt,Hedgehog,Notch" 


# force install decoupler and pertpy to run DE nf
conda activate scanpy-pipeline 
pip install --force-reinstall --no-cache-dir decoupler pertpy
#pip install --force-reinstall --no-cache-dir mira-multiome 

#install
mamba env create -f environment_multiome.yml --override-channels -c conda-forge -c bioconda
mamba env create -f environment_scRNAseq.yml --override-channels -c conda-forge -c bioconda

conda activate scanpy-pipeline
pip install pyscdblfinder

#login node
conda activate multiome-pipeline
# install required packages for magic-impute
pip install scprep --no-deps
pip install decorator
pip install graphtools --no-deps
pip install future
pip install pygsp tasklogger
pip install magic-impute --no-deps
