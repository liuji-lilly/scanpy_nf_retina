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

# 1. READ ANNDATA FIRST
# ==========================================
print("Loading data for ${meta.id}...")
adata = sc.read_h5ad("${adata_annotated}")
# adata = sc.read_h5ad("/fsx/home/l128405/scanpy_nf_retina/results/08_subclustering/test2_adata.3_13_adata_clustered.h5ad")
# # 1. Extract a subset or pointer to the data array
# counts_matrix = adata.layers['counts']

# Retrieving via python
msigdb = dc.op.resource("MSigDB", organism=${organism})')
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


dc.mt.aucell(
    data=adata,
    net=reactome_renamed,
    raw=False,
)

reactome_wnt_pathways = reactome_renamed[
    reactome_renamed["source"].str.contains(${pathway_name}, case=False)
]

# 2. View the unique Wnt pathways you found
wnt_pathways=list(reactome_wnt_pathways["source"].unique())

adata.obs[wnt_pathways] = adata.obsm["score_aucell"][wnt_pathways]

sc.pl.pca(
    adata,
    color=["leiden", "total_counts"] + wnt_pathways,
    frameon=False,
    ncols=2,
    wspace=0.3,
)

fig = plt.gcf()
fig.savefig("wnt_pathway_pca.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

sc.pl.stacked_violin(
    adata,
    var_names=wnt_pathways,
    groupby=["leiden"],
    swap_axes=True,
    figsize=(12, 3),
)

fig = plt.gcf()
fig.savefig("wnt_pathway_stacked_violin.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

progeny = dc.op.progeny(organism=${organism}, top=100)

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

sc.pl.stacked_violin(
    adata,
    var_names=allpathways,
    groupby=["leiden"],
    swap_axes=True,
    figsize=(12, 3),
)

fig = plt.gcf()
fig.savefig("progeny_all_pathways_stacked_violin.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

