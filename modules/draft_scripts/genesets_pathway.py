#!/usr/bin/env python3
"""
Single-Cell Retina Pathway Activation & Cell Type Module Scoring Pipeline
=========================================================================
Authors: Bioinformatics Core
Context: Nextflow Cell Type / Pathway Rendering Module Script Template

Calculates continuous activity scores for:
  1. MSigDB Functional Reactome Pathways (via AUCell)
  2. PROGENy Downstream Signaling Passways (via ULM)
  3. Custom Cell Type Marker Lists (via AUCell)
Outputs optimized publication-quality PDF grids.
"""

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
import math
import csv
# =========================================================================
# 1. DEFINE YOUR PIECEWISE PIPELINE COLOR PALETTE (Resolves NameError)
# =========================================================================
lillycolor = ["#0F3A85","#8A969E","#99BFE5","#E1251B", "#212121", "#FFC709","#511207", "#F58E7D", "#E4EBF1", "#FDE8E5", "#144B2D", "#C6DCD8"]
pal1 = ["#8B633DFF", "#8A999DFF", "#C4B57CFF", "#CB4E33FF", "#F0D1A9FF", "#1B3F8AFF", "#68A1D3FF", "#AECCDFFF", "#D18860FF", "#30588BFF", "#C2944FFF", "#F9D073FF"]
pal2 = ["#C70E7BFF", "#FC6882FF", "#007BC3FF", "#54BCD1FF", "#EF7C12FF", "#F4B95AFF", "#009F3FFF", "#8FDA04FF", "#AF6125FF", "#F4E3C7FF", "#B25D91FF", "#EFC7E6FF", "#EF7C12FF", "#F4B95AFF"]
pal3 = ["#402733FF", "#8F6052FF", "#2F6349FF", "#AE2B47FF", "#7BC17DFF", "#68D8F2FF", "#7982F7FF", "#E03E4AFF", "#FFF759FF", "#ECA987FF", "#705549FF", "#9F4739FF"]
pal4 = ["#3B2319FF", "#80521CFF", "#D29C44FF", "#EBC174FF", "#EDE2CCFF", "#7EC5F4FF", "#4585B7FF", "#225E92FF", "#183571FF", "#43429BFF", "#5E65BEFF"]
pal5 = ["#FF0066FF", "#328C97FF", "#D1AAC2FF", "#A5506DFF", "#B3E0BFFF", "#2A9D3DFF", "#EDF181FF", "#DB7003FF", "#FBA600FF", "#F8C1A6FF", "#A30000FF", "#FF3200FF", "#011A51FF", "#97D1D9FF", "#916C37FF"]

customcolor = lillycolor + pal1 + pal2 + pal3 + pal4 + pal5

# Dynamically add the directory containing sctypepy to the Python Path
sys.path.append(os.path.join("${workflow.projectDir}", "bin"))

sc.settings.verbosity = 3

# 1. READ ANNDATA FIRST
# ==========================================
print("Loading data for ${meta.id}...")
adata = sc.read_h5ad("${adata_annotated}")
# adata = sc.read_h5ad("/fsx/home/l128405/10xgenomics/retina_dataset/Li_2024/06_annotation/retina_P90_adata_annotated.h5ad")
# # 1. Extract a subset or pointer to the data array
# counts_matrix = adata.layers['counts']

# =========================================================================
# msigdb collections download 
# =========================================================================

# Retrieving via python
msigdb = dc.op.resource("MSigDB", organism="${organism}")')
#msigdb = dc.op.resource("MSigDB", organism="mouse")

collections_to_test = {
    "hallmark": "hallmark",
    "kegg": "kegg_pathways",
    "go_bp": "go_biological_process",
    "go_mf": "go_molecular_function",
    "go_cc": "go_cellular_component",
    "pid": "pid_pathways", 
    "wiki": "wikipathways",
    "c3": "transcription_factor_targets",
    "c8": "cell_type_signatures"
}

# Get MSigDB databases
reactome = msigdb.query("collection == 'reactome_pathways'")
go_bp = msigdb.query("collection == 'go_biological_process'")
kegg = msigdb.query("collection == 'kegg_pathways'")
hallmark = msigdb.query("collection == 'hallmark'")

# =========================================================================
# reactome scoring
# =========================================================================

# Filter duplicates
reactome = reactome[~reactome.duplicated(("geneset", "genesymbol"))]

#remove potential noisy gene sets by removing sets with too low (<15) or too high (>500) number of members
geneset_size = reactome.groupby("geneset").size()
reactome = reactome.loc[[15 < geneset_size.loc[gset] < 500 for gset in reactome["geneset"]]]
#reactome

# 1. Rename columns to match what decoupler expects
reactome_renamed = reactome.rename(columns={
    'geneset': 'source', 
    'genesymbol': 'target'
})

# =========================================================================
# STEP 1: GLOBAL REACTOME FUNCTIONAL SCORES (AUCell Engine)
# =========================================================================

dc.mt.aucell(
    data=adata,
    net=reactome_renamed,
    raw=False,
)

# Parse multi-pathway regex inputs passed down via Nextflow core parameters
raw_input = "${params.pathway_names}"
processed_pattern = "|".join([p.strip() for p in raw_input.split(",") if p.strip()])

# Run a regex-safe search against the source column
reactome_selected_pathways = reactome_renamed[
    reactome_renamed["source"].str.contains(processed_pattern, case=False)
]

# View the matching unique pathways found
selected_pathways = list(reactome_selected_pathways["source"].unique())
print(f"Matched the following pathways: {selected_pathways}")

# reactome_selected_pathways = reactome_renamed[
#     reactome_renamed["source"].str.contains("${pathway_name}", case=False)
# ]

# # 2. View the unique Wnt pathways you found
# selected_pathways=list(reactome_selected_pathways["source"].unique())

adata.obs[selected_pathways] = adata.obsm["score_aucell"][selected_pathways]

sc.pl.pca(
    adata,
    color=["leiden", "total_counts"] + selected_pathways,
    frameon=False,
    ncols=2,
    wspace=0.3,
)

fig = plt.gcf()
fig.savefig("selected_pathway_pca.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

num_pathways = len(selected_pathways)
num_clusters = len(adata.obs['leiden'].unique())
custom_width_violin = max(14, num_clusters * 0.5 + 4)
custom_height_violin = max(10, num_pathways * 0.45 + 2)
sc.pl.stacked_violin(
    adata,
    var_names=selected_pathways,
    groupby=["leiden"],
    swap_axes=True,
    figsize=(custom_width_violin, custom_height_violin),
    stripplot=False,              # Keeps the density profile pristine
    show=False,
    width=0.85,  
    density_norm='width',
    dodge=False
)
fig = plt.gcf()
fig.savefig("selected_pathway_stacked_violin.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

# Calculate a clean dynamic canvas size matching your track coordinates
custom_width_dotplot = max(14, num_clusters * 0.5 + 4)
custom_height_dotplot = max(10, num_pathways * 0.45 + 2)

sc.pl.dotplot(
    adata,
    var_names=selected_pathways,
    groupby="leiden",
    dendrogram=True,               # Automatically groups similar clusters by activity profile
    swap_axes=True,                # RETAINED: Keeps clusters on the horizontal X-axis
    expression_cutoff=0.0,         # Hides background noise bubbles below zero activity
    figsize=(custom_width_dotplot, custom_height_dotplot),
    show=False
)

# Export clean high-resolution file
fig = plt.gcf()
fig.savefig("selected_pathway_dotplot.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)
print("Selected pathways dotplot exported successfully.")


# Pull continuous score metrics directly into long-form pandas dataframe
pathway_df = sc.get.obs_df(adata, keys=selected_pathways + ['leiden'])
pathway_df['leiden'] = pd.Categorical(
    pathway_df['leiden'], 
    categories=sorted(pathway_df['leiden'].unique(), key=int), 
    ordered=True
)

# Configure a clean 3-column layout grid matrix
num_paths = len(selected_pathways)
if num_paths <= 2:
    ncols_p = num_paths
elif num_paths == 4:
    ncols_p = 2
else:
    ncols_p = 3  # Standard balanced grid width

nrows_p = math.ceil(num_paths / ncols_p)

panel_width_p = 7.5   
panel_height_p = 4.5   

fig, axes = plt.subplots(
    nrows=nrows_p, 
    ncols=ncols_p, 
    figsize=(ncols_p * panel_width_p, nrows_p * panel_height_p)
)

# Inject spatial horizontal/vertical padding margins between subplots
plt.subplots_adjust(wspace=0.3, hspace=0.5)

if num_paths == 1:
    axes = [axes]
else:
    axes = axes.flatten()

# Slice colors to match the exact number of active clusters
cluster_colors = customcolor[:len(pathway_df['leiden'].cat.categories)]

# Loop through each pathway and generate an un-stacked widescreen layout
for i, pathway in enumerate(selected_pathways):
    sns.violinplot(
        data=pathway_df,
        x='leiden',
        y=pathway,
        ax=axes[i],
        inner=None,             # Removes internal box markers to keep shapes clean
        density_norm='width',    # Standardizes distribution heights across clusters
        hue='leiden',           # Colors distributions by cluster
        palette=cluster_colors, # Synchronizes colors with your global pipeline
        dodge=False,            # CRITICAL: Snaps violins perfectly onto tick centers
        width=0.85,             # Standardized optimal track width per lane
        legend=False
    )
    # Clean title text formatting (Strips away raw database prefixes for clean titles)
    clean_title = pathway.replace("REACTOME_", "").replace("_", " ")
    if len(clean_title) > 40:
        clean_title = clean_title[:40] + "..." # Protects layouts from ultra-long text block
    axes[i].set_title(clean_title, fontsize=12, fontweight='bold', pad=12)
    axes[i].set_xlabel('Leiden Clusters', fontsize=11, labelpad=6)
    axes[i].set_ylabel('AUCell Activity Score', fontsize=10)
    # Keep cluster numbers flat, horizontal, and perfectly aligned under points
    axes[i].tick_params(axis='x', labelsize=10, rotation=0)
    axes[i].set_xticklabels(axes[i].get_xticklabels(), ha='center')
    
    # Add light background layout grids for easier score mapping
    axes[i].grid(True, axis='y', linestyle='--', alpha=0.3)

# Cleanly remove empty plot backdrops from the trailing edge of the grid matrix
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Save vector high-resolution clear file
plt.savefig("selected_pathway_individual_violins.pdf", bbox_inches='tight', dpi=300)
plt.close()
print("Individual un-stacked pathway violin panel views exported successfully.")


# =========================================================================
# Progeny Pathway Scoring 
# =========================================================================

progeny = dc.op.progeny(organism="${organism}", top=100)
#progeny = dc.op.progeny(organism="mouse", top=100)

dc.mt.ulm(
    data=adata,
    net=progeny,
    raw=False,
)

allpathways=adata.obsm["score_ulm"].columns.tolist()
adata.obs[allpathways] = adata.obsm["score_ulm"][allpathways]

sc.pl.pca(
    adata,
    color=["leiden", "total_counts"] + allpathways,
    frameon=False,
    ncols=2,
    wspace=0.3,
)

fig = plt.gcf()
fig.savefig("progeny_all_pathways_pca.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

num_pathways = len(allpathways)
num_clusters = len(adata.obs['leiden'].unique())
custom_width_violin = max(14, num_clusters * 0.5 + 4)
custom_height_violin = max(10, num_pathways * 0.45 + 2)

sc.pl.stacked_violin(
    adata,
    var_names=allpathways,
    groupby=["leiden"],
    swap_axes=True,
    figsize=(custom_width_violin, custom_height_violin),
    stripplot=False,              # Keeps the density profile pristine
    show=False,
    density_norm='width',
    dodge=False                  # <--- CRITICAL FIX: Forces the violin envelope to expand horizontally
)

fig = plt.gcf()
fig.savefig("progeny_all_pathways_stacked_violin.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

custom_width_dotplot = max(12, len(allpathways) * 0.5 + 4)
custom_height_dotplot = max(8, len(adata.obs['leiden'].unique()) * 0.45 + 2)

sc.pl.dotplot(
    adata,
    var_names=allpathways,
    groupby="leiden",
    dendrogram=True,               # Hierarchical clustering of rows based on pathway activity
    swap_axes=False,               # Keeps pathways on rows, clusters on columns to match standard layouts
    expression_cutoff=0.0,         # Minimum value threshold for background dot suppression
    figsize=(custom_width_dotplot, custom_height_dotplot),
    show=False
)

# Export clean high-resolution file
fig = plt.gcf()
fig.savefig("progeny_all_pathways_dotplot.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

# =========================================================================
# custom marker gene network scoring
# =========================================================================

#marker_genes_file="/fsx/home/l128405/10xgenomics/retina_dataset/retina_markergenes.csv"
network_rows = []

try:
    with open("${celltyp_markergenes_file}", mode="r", encoding="utf-8-sig") as file:
        file.readline()  # Skip header row if it exists
        file.seek(0)
        reader = csv.reader(file)
        rows = list(reader)
except UnicodeDecodeError:
    with open("${celltyp_markergenes_file}", mode="r", encoding="cp1252") as file:
        reader = csv.reader(file)
        rows = list(reader)

# Convert horizontal CSV rows into long-form pairwise entries
for row in rows:
    if row:
        cell_type = row[0].strip()
        genes = [g.strip() for g in row[1:] if g.strip()]
        
        # Verify gene presence in your specific dataset to prevent AUCell runtime errors
        valid_genes = [g for g in genes if g in adata.var_names]
        
        for gene in valid_genes:
            network_rows.append({'source': cell_type, 'target': gene})

# Create the final network DataFrame required by Decoupler
marker_network = pd.DataFrame(network_rows)
print(f"Successfully built network with {len(marker_network['source'].unique())} cell types and {len(marker_network)} total interactions.")

# =========================================================================
# 2. RUN AUCELL OVER THE CELL TYPE MODULES
# =========================================================================
print("Running AUCell module scoring for all cell types...")
dc.mt.aucell(
    data=adata,
    net=marker_network,
    raw=False,
    tmin=5 # Minimum number of genes required to compute AUCell score
)
# Extract the unique cell type names to handle plotting arrays
#cell_type_modules = list(marker_network["source"].unique())
cell_type_modules = list(adata.obsm["score_aucell"].columns)

# Map the matrix results back into main adata.obs for easy Scanpy access
adata.obs[cell_type_modules] = adata.obsm["score_aucell"][cell_type_modules]


# Plot all cell-type AUCell scores in a balanced 3-column unstacked violin grid
num_modules = len(cell_type_modules)
ncols = 3
nrows = math.ceil(num_modules / ncols)

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 7.0, nrows * 4.0))
plt.subplots_adjust(wspace=0.3, hspace=0.45)
axes = axes.flatten()

for i, cell_type in enumerate(cell_type_modules):
    sns.violinplot(
        data=sc.get.obs_df(adata, keys=[cell_type, 'leiden']),
        x='leiden',
        y=cell_type,
        ax=axes[i],
        inner=None,
        density_norm='width',
        hue='leiden',
        dodge=False, # Keeps violins centered on x-ticks
        width=0.85,
        legend=False
    )
    axes[i].set_title(f"{cell_type} (AUCell)", fontsize=14, fontweight='bold')
    axes[i].set_xlabel('Leiden Clusters', fontsize=11)
    axes[i].set_ylabel('Activity Score', fontsize=10)
    axes[i].tick_params(axis='x', labelsize=10)

# Hide any trailing blank frames
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.savefig("celltype_aucell_individual_violins.pdf", bbox_inches='tight')
plt.close()


# =========================================================================
# 2. PREPARE DATAFRAME COLUMNS & SANITIZE KEYS
# =========================================================================
# Extract computed cell type columns
cell_type_modules = list(adata.obsm["score_aucell"].columns)

# Map them back into your observations metadata matrix slot safely
if cell_type_modules:
    adata.obs[cell_type_modules] = adata.obsm["score_aucell"][cell_type_modules]
else:
    raise ValueError("No cell types passed your tmin=3 filtering threshold!")

# Dynamically pull the exact keys that exist to prevent KeyError
swarm_df = sc.get.obs_df(adata, keys=cell_type_modules + ['leiden']).copy()

# Force standard string representation on clusters
swarm_df['leiden'] = swarm_df['leiden'].astype(str)

# Ensure cell type headers contain no spaces or illegal characters for dataframe interpretation
# This translates 'Amacrine Cells' into 'Amacrine_Cells' safely
rename_map = {col: col.replace(" ", "_").replace("-", "_") for col in cell_type_modules}
swarm_df = swarm_df.rename(columns=rename_map)
sanitized_modules = [rename_map[col] for col in cell_type_modules]

# Sort cluster names numerically
unique_clusters = sorted(swarm_df['leiden'].unique(), key=lambda x: int(x) if x.isdigit() else 999)
cluster_colors = customcolor[:len(unique_clusters)]

# =========================================================================
# 3. CONFIGURE MULTI-PANEL GRID
# =========================================================================
num_panels = len(sanitized_modules)
ncols_s = 3
nrows_s = math.ceil(num_panels / ncols_s)

panel_width_s = 7.5   
panel_height_s = 4.5   

fig, axes = plt.subplots(
    nrows=nrows_s, 
    ncols=ncols_s, 
    figsize=(ncols_s * panel_width_s, nrows_s * panel_height_s)
)

plt.subplots_adjust(wspace=0.3, hspace=0.5)

if num_panels == 1:
    axes = [axes]
else:
    axes = axes.flatten()

# Track active panel assignments dynamically to avoid gaps from skipped slots
active_plots = 0

# =========================================================================
# 4. PLOT TRACK DISTRIBUTION LOOPS
# =========================================================================
for cell_type in sanitized_modules:
    # CRITICAL SAFETY GUARD: Skip if the module name is missing from the data frame columns
    if cell_type not in swarm_df.columns:
        print(f"Skipping plot for '{cell_type}': column not found in data frame slice.")
        continue
        
    sns.stripplot(
        data=swarm_df,
        x='leiden',
        y=cell_type,
        ax=axes[active_plots],
        hue='leiden',           # Colors points by cluster identity
        palette=cluster_colors, # Synchronizes colors across your pipeline
        order=unique_clusters,  # Forces strict numerical sorting on x-axis
        hue_order=unique_clusters,
        jitter=0.25,            # Spreads points horizontally to show density
        size=2.5,               # Crisp point diameter prevents overlap saturation
        dodge=False,            # Centered directly over the numerical tick line
        alpha=0.65,             # Layering transparency reveals true hub density
        legend=False
    )
    
    # Subplot text properties
    clean_title = cell_type.replace('_', ' ')
    axes[active_plots].set_title(f"{clean_title} Profile", fontsize=13, fontweight='bold', pad=12)
    axes[active_plots].set_xlabel('Leiden Clusters', fontsize=11, labelpad=6)
    axes[active_plots].set_ylabel('AUCell Activity Score', fontsize=10)
    
    # Keep tick elements horizontal and centered perfectly under point tracks
    axes[active_plots].tick_params(axis='x', labelsize=10, rotation=0)
    axes[active_plots].set_xticklabels(axes[active_plots].get_xticklabels(), ha='center')
    axes[active_plots].grid(True, axis='y', linestyle='--', alpha=0.3)
    
    active_plots += 1

# Hide any unused trailing subplots cleanly using the visibility flag
for j in range(active_plots, len(axes)):
    axes[j].set_visible(False)

# =========================================================================
# 5. RENDER AND EXPORT CLEAN DOCUMENT
# =========================================================================
plt.savefig("celltype_aucell_individual_stripplots.pdf", bbox_inches="tight", dpi=300)
plt.close()
print(f"Point distribution grids exported successfully for {active_plots} cell types.")