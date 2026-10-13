#!/usr/bin/env python3
import sys
import os
import glob
import argparse
import scanpy as sc
import pandas as pd
from datetime import datetime

def main():
    # Set up robust command-line argument parsing
    parser = argparse.ArgumentParser(description="Generate Master Summary Report for scRNA-seq Pipeline")
    parser.add_argument("--adata_dir", required=True, help="Directory containing annotated h5ad files")
    parser.add_argument("--qc_dir", required=True, help="Directory containing QC metrics CSV files")
    parser.add_argument("--marker_dir", required=True, help="Directory containing marker gene CSV files")
    args = parser.parse_args()

    # Parse directories using the verified arguments
    adata_files = glob.glob(os.path.join(args.adata_dir, "*_adata_annotated.h5ad"))
    qc_files = glob.glob(os.path.join(args.qc_dir, "*_qc_metrics.csv"))
    markers_files = glob.glob(os.path.join(args.marker_dir, "*.csv"))

    # Map files by cleanly pulling out the core sample prefix
    adata_map = {os.path.basename(f).replace("_adata_annotated.h5ad", ""): f for f in adata_files}
    qc_map = {os.path.basename(f).replace("_qc_metrics.csv", ""): f for f in qc_files}
    markers_map = {os.path.basename(f).split('_')[0]: f for f in markers_files}

    # Identify shared matching sample IDs across components
    all_samples = sorted(list(set(adata_map.keys()) & set(qc_map.keys())))

    if not all_samples:
        print("Error: No matching sample IDs found between h5ad and QC metrics files.")
        print(f"Looked inside adata_dir: '{args.adata_dir}' -> Found files: {adata_files}")
        print(f"Looked inside qc_dir: '{args.qc_dir}' -> Found files: {qc_files}")
        print(f"Parsed AnnData IDs: {list(adata_map.keys())}")
        print(f"Parsed QC Metric IDs: {list(qc_map.keys())}")
        sys.exit(1)

    # Initialize Master Report Content
    md = f"# Multi-Sample Single-Cell RNA-seq Summary Report\n\n"
    md += f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
    md += f"**Total Samples Processed:** {len(all_samples)}\n\n"

    # --- SECTION 1: MASTER OVERVIEW TABLE ---
    md += "## Master Dataset Statistics\n\n"
    md += "| Sample ID | Total Cells | Total Genes | Clusters | Batch | Treatment |\n"
    md += "|---|---|---|---|---|---|\n"

    sample_details = {} 

    for sample_id in all_samples:
        adata = sc.read_h5ad(adata_map[sample_id])
        n_cells = adata.n_obs
        n_genes = adata.n_vars
        n_clusters = adata.obs['leiden'].nunique() if 'leiden' in adata.obs.columns else 0
        batch = adata.obs['batch'].iloc[0] if 'batch' in adata.obs.columns else "N/A"
        treatment = adata.obs['treatment'].iloc[0] if 'treatment' in adata.obs.columns else "N/A"

        md += f"| **{sample_id}** | {n_cells:,} | {n_genes:,} | {n_clusters} | {batch} | {treatment} |\n"
        
        sample_details[sample_id] = {
            'adata': adata, 
            'n_cells': n_cells, 
            'n_clusters': n_clusters
        }

    md += "\n"

    # --- SECTION 2: INDIVIDUAL SAMPLE BREAKDOWNS ---
    for sample_id in all_samples:
        md += f"--- \n\n## Sample Breakdowns: {sample_id}\n\n"
        
        adata = sample_details[sample_id]['adata']
        n_cells = sample_details[sample_id]['n_cells']
        
        # QC Metrics Breakdown
        md += "### QC Metrics\n\n"
        qc_df = pd.read_csv(qc_map[sample_id])
        
        # Normalize casing
        qc_df.columns = [str(c).lower() for c in qc_df.columns]
        
        metric_col = 'metric' if 'metric' in qc_df.columns else qc_df.columns[0]
        val_col = 'median' if 'median' in qc_df.columns else ('value' if 'value' in qc_df.columns else qc_df.columns[1])
        
        md += f"| Metric | Summary Value (1 decimal place) |\n|---|---|\n"
        for _, row in qc_df.iterrows():
            try:
                raw_val = float(row[val_col])
                val_str = f"{raw_val:,.1f}"
            except (ValueError, TypeError, KeyError):
                val_str = str(row.get(val_col, 'N/A'))
            md += f"| {row[metric_col]} | {val_str} |\n"
        md += "\n"

        # Cell Type / Cluster Distribution
        if 'sctype_classification' in adata.obs.columns or 'cell_type' in adata.obs.columns:
            col = 'sctype_classification' if 'sctype_classification' in adata.obs.columns else 'cell_type'
            md += "### Cell Type Distribution\n\n| Cell Type | Cell Count | Percentage |\n|---|---|---|\n"
            for ct, count in adata.obs[col].value_counts().items():
                pct = (count / n_cells) * 100
                md += f"| {ct} | {count:,} | {pct:.1f}% |\n"
        else:
            md += "### Cluster Distribution\n\n| Cluster | Cell Count |\n|---|---|\n"
            if 'leiden' in adata.obs.columns:
                for cl, count in adata.obs['leiden'].value_counts().sort_index().items():
                    md += f"| {cl} | {count:,} |\n"
        md += "\n"

        # Top Marker Genes
        if sample_id in markers_map:
            md += "### Top Marker Genes per Cluster\n\n"
            markers_df = pd.read_csv(markers_map[sample_id])
            if 'cluster' in markers_df.columns:
                gene_col = 'names' if 'names' in markers_df.columns else markers_df.columns[0]
                for cluster in markers_df['cluster'].unique():
                    cluster_markers = markers_df[markers_df['cluster'] == cluster].head(5)
                    md += f"**Cluster {cluster}:** " + ", ".join(cluster_markers[gene_col].astype(str).tolist()) + "\n\n"

    # --- SAVE OUTPUTS ---
    with open("combined_analysis_summary.md", "w") as f:
        f.write(md)

    html = f"<!DOCTYPE html><html><head><title>Multi-Sample Analysis Report</title>" \
           f"<style>body{{font-family:Arial,sans-serif;margin:40px;line-height:1.6;color:#333;}}" \
           f"h1, h2{{color:#2c3e50;border-bottom:2px solid #3498db;padding-bottom:5px;}}h3{{color:#2980b9;}}" \
           f"table{{border-collapse:collapse;width:100%;margin-bottom:30px;box-shadow:0 2px 3px rgba(0,0,0,0.1);}}" \
           f"th,td{{border:1px solid #ddd;padding:12px;text-align:left;}}th{{background:#3498db;color:#fff;}}tr:nth-child(even){{background:#f9f9f9;}}</style>" \
           f"</head><body>{md}</body></html>"

    with open("combined_analysis_summary.html", "w") as f:
        f.write(html)
    print("Reports generated successfully.")

if __name__ == "__main__":
    main()