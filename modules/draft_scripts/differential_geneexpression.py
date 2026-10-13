#!/usr/bin/env python3
import decoupler as dc
import pertpy as pt
import lamindb as ln
import numpy as np
import pertpy as pt
import scanpy as sc
import pandas as pd
import anndata
import re
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

# Dynamically add the directory containing sctypepy to the Python Path
sys.path.append(os.path.join("${workflow.projectDir}", "bin"))

sc.settings.verbosity = 3

import run_pathway_enrichment
from run_pathway_enrichment import run_pathway_enrichment


# 1. READ ANNDATA FIRST
# ==========================================
print("Loading data for ${meta.id}...")
adata = sc.read_h5ad("${adata_annotated}")
# adata = sc.read_h5ad("/fsx/home/l128405/scanpy_nf_retina/results/08_subclustering/test2_adata.3_13_adata_clustered.h5ad")
# # 1. Extract a subset or pointer to the data array
# counts_matrix = adata.layers['counts']

# # 2. Extract raw array data (handles dense numpy arrays or sparse scipy matrices)
# if hasattr(counts_matrix, "data"):
#     X_counts = counts_matrix.data
# else:
#     X_counts = counts_matrix

# # 3. Check if all values equal their rounded integers
# is_integer = np.array_equal(X_counts, np.round(X_counts))
# print(f"Is the counts layer strictly integers? {is_integer}")

# # filter cells which have less than 200 genes 
# # filter genes which were found in less than 3 cells for a rudimentary quality control
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)

# ==========================================
# 2. READ CSV & EXTRACT CELLS OF INTEREST
# ==========================================
print("Reading target cells of interest from CSV...")
df = pd.read_csv('${DE_csv}')
#df = pd.read_csv('/fsx/home/l128405/scanpy_nf_retina/bin/DE_samplesheet_pseudobulk.csv')

df['target_identities'] = df['target_identities'].astype(str)

adata_subsets = []

for h5ad_path, group_df in df.groupby('input_h5ad'):
    print(f"\nLoading file: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    
    for idx, row in group_df.iterrows():
        condition = row["condition"]
        col_info = row["col_info"]
        target_id_raw = row["target_identities"]
        
        # Clean up rep input string values
        rep_val = str(row["rep"]).strip()
        rep_col_name = 'NA' if rep_val.upper() in ['NA', 'NAN', ''] else rep_val
        
        target_ids_list = [item.strip() for item in target_id_raw.split(',')]
        print(f" -> Extracting cells for '{condition}' ({col_info} matching {target_ids_list})...")
        
        adata_subset = adata[adata.obs[col_info].astype(str).isin(target_ids_list)].copy()
        
        if adata_subset.n_obs == 0:
            print(f"    ⚠️ Warning: Zero cells found matching {col_info} in {target_ids_list}!")
            continue
            
        # Unified tracking attributes
        adata_subset.obs["test_group"] = condition.lower() 
        adata_subset.obs["rep_col_name"] = rep_col_name
        
        # Build replicate group identifiers only if a valid column is provided
        if rep_col_name != 'NA':
            adata_subset.obs["rep_group"] = adata_subset.obs[rep_col_name].astype(str) + "_" + condition.lower()
        else:
            adata_subset.obs["rep_group"] = 'NA'
            
        adata_subsets.append(adata_subset)

if adata_subsets:
    # Concatenate and verify
    adata_combined = sc.concat(adata_subsets, join="outer")
    print("\n--- Summary of combined dataset ---")
    print(adata_combined)


# 3. Check replicate configuration state to route execution path
# Look at unique rep entries (ignoring 'NA' designations)
unique_reps = [r for r in adata_combined.obs["rep_col_name"].unique() if r != 'NA']

if len(unique_reps) == 0:
    # =========================================================================
    # TRACK A: CELL LEVEL SINGLE-CELL DEG (Standard Wilcoxon Engine)
    # =========================================================================
    print("\n[Track: Single-Cell] No biological replicate tracking requested. Executing cell-level Wilcoxon...")
    # FIX: Create normalized copy ONLY for the single-cell test track
    adata_deg = adata_combined.copy()
    # 1. Filter structural outliers cleanly
    #sc.pp.filter_cells(adata_deg, min_counts=1)
    #sc.pp.filter_genes(adata_deg, min_cells=1)
    adata_deg.X = adata_deg.layers["counts"].copy()
    # 2. Normalize and log-transform the data (Critical to prevent NaNs!)
    sc.pp.normalize_total(adata_deg, target_sum=1e4)
    sc.pp.log1p(adata_deg)

    # 2. Run Wilcoxon directly on the existing combined layer
    sc.tl.rank_genes_groups(
        adata_deg, 
        groupby='test_group', 
        groups=['test'], 
        reference='control', 
        method='wilcoxon',
        use_raw=False
    )

    # Extract the target dataframe
    deg_df = sc.get.rank_genes_groups_df(adata_deg, group='test')
    deg_df.to_csv('single_cell_test_vs_control_deg.csv', index=False)
    pds2_df = deg_df.rename(
        columns={
            "names": "variable",
            "logfoldchanges": "log_fc",
            "pvals_adj": "adj_p_value",
        }
    )
    # run pathway enrichment and save outputs
    run_pathway_enrichment(deg_dataframe=pds2_df, prefix="singlecell", organism=${organism})

    # Apply your significance filter and plot safely
    filtered_df = pds2_df[pds2_df["adj_p_value"] < 0.01].copy()
    pds2.plot_fold_change(filtered_df, n_top_vars=15)
    # Save the figure
    fig = plt.gcf()
    fig.savefig("singlecell_DE_barplot.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    # Generate volcano plot setting fold-change line threshold at 0
    pds2.plot_volcano(pds2_df, log2fc_threshold=0)
    # Target active canvas and export
    fig = plt.gcf()
    fig.savefig("singlecell_DE_volcano.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print("Single-cell DEG execution completed successfully.")

    

else:
    # =========================================================================
    # TRACK B: PSEUDOBULK SAMPLING DEG
    # =========================================================================
    rep_column = unique_reps[0] # Dynamic recovery of your true sample column ID
    print(f"\n[Track: Pseudobulk] Replicate metadata identified: '{rep_column}'. Running Pseudobulk...")
        # Aggregate cell expression profiles up to individual biological sample nodes

    # 1. Create unique pseudobulk group labels
    adata_combined.obs["sample"] = adata_combined.obs["rep_group"]
    # 2. pseudobulk aggregation using decoupler
    adata_pb = dc.pp.pseudobulk(
        adata_combined,
        sample_col="sample",
        groups_col=None, 
        layer='counts',
        mode="sum"
    )

    # 3. Filter low-quality samples (Keep samples with >=10 cells and >=1000 counts)
    dc.pp.filter_samples(adata_pb, min_cells=10,min_counts=1000)

    # 4. Save raw counts layer for DESeq2 before transformations
    adata_pb.layers['counts'] = adata_pb.X.copy()

    # 5. Transform data for PCA visualization tracking
    sc.pp.normalize_total(adata_pb, target_sum=1e6)
    sc.pp.log1p(adata_pb)
    sc.pp.scale(adata_pb, max_value=10)
    sc.tl.pca(adata_pb)

    # 6. Swap layer back to raw counts (Critical for DESeq2!)
    dc.pp.swap_layer(adata=adata_pb, key="counts", inplace=True)

    # 7. Add log metrics for metadata tracking
    adata_pb.obs["psbulk_counts_log"] = np.log(adata_pb.obs["psbulk_counts"])
    adata_pb.obs["psbulk_cells_log"] = np.log(adata_pb.obs["psbulk_cells"])

    # 8. Run PCA association and save plots
    valid_keys = [col for col in adata_pb.obs.columns if adata_pb.obs[col].nunique() > 1]
    dc.tl.rankby_obsm(adata_pb, key="X_pca", obs_keys=valid_keys)
    sc.pl.pca_variance_ratio(adata_pb, log=True, n_pcs=20, show=False)
    plt.savefig('pca_variance_ratio.pdf', bbox_inches='tight')
    plt.close()

    # 9. Save out combined PCA breakdown plot
    sc.pl.pca(adata_pb, color=sorted(list(adata_pb.obs.columns)), ncols=2, size=300, show=False)
    fig = plt.gcf()
    fig.savefig("pseudobulk_pca.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    # 10. Filter out low-expression genes before running PyDESeq2
    dc.pp.filter_by_expr(
        adata=adata_pb,
        group="test_group",
        min_count=10,
        min_total_count=15,
        large_n=10,
        min_prop=0.7,
    )
    dc.pp.filter_by_prop(
        adata=adata_pb,
        min_prop=0.1,
        min_smpls=2,
    )


    # 11. Initialize PyDESeq2 on the perfectly cleaned matrix
    pds2 = pt.tl.PyDESeq2(adata=adata_pb, design='~test_group')
    pds2.fit()

    res_df = pds2.test_contrasts(pds2.contrast(
        column="test_group",
        baseline="control",
        group_to_compare="test")
    )
    res_df.to_csv("pseudobulk_DE_results.csv", index=True)

    #run pathway enrichment and save outputs
    run_pathway_enrichment(deg_dataframe=res_df, prefix="pseudobulk", organism=${organism})

    pds2.plot_volcano(res_df, log2fc_threshold=0)
    fig = plt.gcf()
    fig.savefig("pseudobulk_DE_volcano.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)


    pds2.plot_fold_change(res_df[res_df["adj_p_value"] < 0.01].copy(), n_top_vars=15)
    fig = plt.gcf()
    fig.savefig("pseudobulk_DE_barplot.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    pds2.plot_paired(
        adata_pb,
        results_df=res_df,
        n_top_vars=10,
        groupby="test_group",
        pairedby="sctype_classification"
    )

    fig = plt.gcf()
    fig.savefig("pseudobulk_DE_paired_top10.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)


