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

nextflow run visualization_main.nf \
--custom_marker_file /fsx/home/l128405/scanpy_nf_retina/sampledata/marker_genes.csv   \
--annotate_celltypes false   --outdir "./results"

nextflow run subclustering_main.nf \
  -c nextflow.config \
  --interest_csv "/fsx/home/l128405/scanpy_nf_retina/bin/subclustering.csv" \
  --sub_input_dir "/fsx/home/l128405/scanpy_nf_retina/sampledata/results/06_annotation/" \
  --custom_marker_file "/fsx/home/l128405/scanpy_nf_retina/sampledata/marker_genes.csv" \
  --tissue_name "Mouse Organoid"

# force install decoupler and pertpy to run DE nf
conda activate scanpy-pipeline 
pip install --force-reinstall --no-cache-dir decoupler pertpy
