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

def run_pathway_enrichment(deg_dataframe, prefix="track_output", organism="mouse"):
    """
    Runs automated GSEA (MSigDB), signaling footprint (PROGENY), and 
    Transcription Factor enrichment (CollecTRI) analysis workflows.
    Includes leading-edge pathway heatmaps and TF-to-Pathway footprint link analysis.
    """

    print(f"\n[{prefix.upper()}] Starting Pathway Enrichment & Pathway Inference Processing ({organism})...")
    # 1. Clean the input dataframe and ensure gene symbols are in the index
    res_df = deg_dataframe.copy()
    if "variable" in res_df.columns:
        res_df = res_df.set_index("variable")

    # 2. Dynamically extract structural columns
    stat_col = next((c for c in ["stat", "log_fc", "logFC", "log2FoldChange"] if c in res_df.columns), None)
    padj_col = next((c for c in ["padj", "p_val_adj", "adj.P.Val", "fdr", "pvalue", "p_val", "adj_p_value"] if c in res_df.columns), None)
    if not stat_col or not padj_col:
        print(f"Skipping pathway analysis: Required metrics missing from columns {list(res_df.columns)}")
        return
    
    # 3. Create dynamic contrast label based on output file prefix
    contrast_name = f"{prefix}_contrast"
    data = res_df[[stat_col]].T.rename(index={stat_col: contrast_name})
    pathway_summary_records = []

    # 4. Setup MSigDB pathway resource dynamically by organism variable
    msigdb = dc.op.resource("MSigDB", organism=organism)
    collections_to_test = {
        "hallmark": "hallmark", "kegg": "kegg_pathways", "reactome": "reactome_pathways",
        "go_bp": "go_biological_process", "go_mf": "go_molecular_function", "go_cc": "go_cellular_component"
    }
    #Track variables to systematically isolate the single most significant pathway across all databases for leading-edge plotting
    top_pathway_name, top_pathway_p, top_pathway_genes = None, 1.0, []

    # 5. Loop across MSigDB database paths
    for label, db_name in collections_to_test.items():
        net_df = msigdb.query(f"collection == '{db_name}'")
        # For mouse data, convert uppercase MSigDB human gene symbols to standard Title Case (e.g., SOX9 to Sox9) to enable key matching
        if organism.lower() == "mouse":
            net_df["genesymbol"] = net_df["genesymbol"].str.capitalize()
        net_df = net_df[~net_df.duplicated(("geneset", "genesymbol"))]
        net_renamed = net_df.rename(columns={"geneset": "source", "genesymbol": "target"})
        try:
            # Execute Gene Set Enrichment Analysis (GSEA) statistical engine
            score, padj = dc.mt.gsea(data=data, net=net_renamed)
        except AssertionError:
            continue
        significant_mask = (padj.T < 0.05).iloc[:, 0]
        significant_scores = score.loc[:, significant_mask]
        if not significant_scores.empty:
            fig = dc.pl.barplot(data=significant_scores, name=contrast_name, top=25, return_fig=True)
            fig.subplots_adjust(left=0.45)
            plt.title(f"Top Enriched Pathways - {label.upper()}")
            # FIX: Change bbox_inches to None to prevent triggering internal tight layout warnings
            plt.savefig(f"{prefix}_gsea_barplot_{label}.pdf", bbox_inches=None, format="pdf")
            plt.close(fig)
        available_genes_net = net_renamed[net_renamed["target"].isin(res_df.index)]
        for pathway in score.columns:
            p_adj = padj.at[contrast_name, pathway]
            pathway_genes = available_genes_net[available_genes_net["source"] == pathway]["target"]
            gene_stats = res_df.loc[pathway_genes]
            significant_de_genes = gene_stats[gene_stats[padj_col] < 0.05]
            up_signif = list(significant_de_genes[significant_de_genes[stat_col] > 0].index)
            down_signif = list(significant_de_genes[significant_de_genes[stat_col] < 0].index)
            if p_adj < top_pathway_p and len(significant_de_genes) > 2:
                top_pathway_p = p_adj
                top_pathway_name = pathway
                top_pathway_genes = list(significant_de_genes.index)
            pathway_summary_records.append({
                "Database_Collection": label.upper(), "Pathway_Name": pathway, "Enrichment_Score": score.at[contrast_name, pathway],
                "Pathway_Adj_P_Value": p_adj, "Pathway_Significance": "Significant" if p_adj < 0.05 else "Not Significant",
                "Upregulated_Signif_Genes": ", ".join(up_signif), "Downregulated_Signif_Genes": ", ".join(down_signif),
                "Count_Upregulated_Signif": len(up_signif), "Count_Downregulated_Signif": len(down_signif), "Total_Pathway_Genes_In_Data": len(gene_stats),
            })
    if pathway_summary_records:
        pd.DataFrame(pathway_summary_records).sort_values(by="Pathway_Adj_P_Value", ascending=True).to_csv(f"{prefix}_gsea_summary.csv", index=False)
    if top_pathway_name and top_pathway_genes:
        print(f"[{prefix.upper()}] Generating leading-edge heatmap for top pathway: {top_pathway_name}")
        heatmap_data = res_df.loc[top_pathway_genes, [stat_col]].sort_values(by=stat_col, ascending=False)
        import seaborn as sns
        plt.figure(figsize=(6, max(4, len(top_pathway_genes) * 0.25)))
        sns.heatmap(heatmap_data, cmap="vlag", center=0, annot=True, fmt=".2f", cbar_kws={'label': stat_col})
        plt.title(f"Leading-Edge Genes\n{top_pathway_name}")
        plt.savefig(f"{prefix}_top_pathway_leading_edge_heatmap.pdf", bbox_inches=None)
        plt.close()
    
    print(f"[{prefix.upper()}] Extracting PROGENY signaling topologies...")
    progeny = dc.op.progeny(organism=organism, top=500)
    score_prog, padj_prog = dc.mt.ulm(data=data, net=progeny)
    if not score_prog.loc[:, (padj_prog.T < 0.05).iloc[:, 0]].empty:
        fig = dc.pl.barplot(data=score_prog.loc[:, (padj_prog.T < 0.05).iloc[:, 0]], name=contrast_name, return_fig=True)
        fig.subplots_adjust(left=0.35)
        fig.savefig(f"{prefix}_progeny_barplot.pdf", bbox_inches=None, dpi=300)
        plt.close(fig)
    
    print(f"[{prefix.upper()}] Extracting CollecTRI Transcription Factor regulons...")
    # FIX: Correct module endpoint to map your specific version environment layout
    collectri = dc.op.collectri(organism=organism)
    try:
        score_tf, padj_tf = dc.mt.ulm(data=data, net=collectri)
        if not score_tf.loc[:, (padj_tf.T < 0.05).iloc[:, 0]].empty:
            fig = dc.pl.barplot(data=score_tf.loc[:, (padj_tf.T < 0.05).iloc[:, 0]], name=contrast_name, top=25, return_fig=True)
            fig.subplots_adjust(left=0.35)
            plt.title(f"Top Active Transcription Factors - CollecTRI")
            fig.savefig(f"{prefix}_collectri_tf_barplot.pdf", bbox_inches=None, dpi=300)
            plt.close(fig)
            pd.DataFrame({"TF": score_tf.columns, "Activity_Score": score_tf.loc[contrast_name].values, "Adj_P_Value": padj_tf.loc[contrast_name].values}).sort_values(by="Adj_P_Value", ascending=True).to_csv(f"{prefix}_collectri_tf_summary.csv", index=False)
        sig_tfs = padj_tf.columns[(padj_tf.loc[contrast_name] < 0.05)]
        if len(sig_tfs) > 0:
            print(f"[{prefix.upper()}] Exporting {len(sig_tfs)} significant transcription factors to CSV...")
            # Extract statistics specifically for the significant TFs
            tf_export_df = pd.DataFrame({
                "TF_Name": sig_tfs,
                "Activity_Score": score_tf.loc[contrast_name, sig_tfs].values,
                "Adj_P_Value": padj_tf.loc[contrast_name, sig_tfs].values
            })
            # Sort by lowest adjusted p-value so the most significant TFs are at the top
            tf_export_df = tf_export_df.sort_values(by="Adj_P_Value", ascending=True)
            # Export to a structured CSV tracked by your execution track prefix
            tf_export_df.to_csv(f"{prefix}_significant_transcription_factors.csv", index=False)
        sig_pathways = padj_prog.columns[(padj_prog.loc[contrast_name] < 0.05)]
        if len(sig_tfs) > 0 and len(sig_pathways) > 0:
            print(f"[{prefix.upper()}] Cross-linking active TFs to PROGENY footprint downstream nodes...")
            overlap_records = []
            for tf in sig_tfs[:10]:
                tf_targets = set(collectri[collectri["source"] == tf]["target"])
                for path in sig_pathways:
                    path_targets = set(progeny[progeny["source"] == path]["target"])
                    shared = tf_targets.intersection(path_targets)
                    if len(shared) > 0:
                        overlap_records.append({"TF": tf, "PROGENY_Pathway": path, "Shared_Target_Count": len(shared)})
            if overlap_records:
                link_df = pd.DataFrame(overlap_records).pivot(index="TF", columns="PROGENY_Pathway", values="Shared_Target_Count").fillna(0)
                plt.figure(figsize=(8, 6))
                import seaborn as sns
                sns.heatmap(link_df, cmap="Purples", annot=True, fmt="g")
                plt.title("Upstream TFs to Downstream Pathway Footprint Intersections")
                plt.savefig(f"{prefix}_tf_to_progeny_links.pdf", bbox_inches=None)
                plt.close()
    except AssertionError:
        print(f"Skipping network integrations: Feature vector mismatch.")
    print(f"[{prefix.upper()}] Workflow pipeline completed execution successfully.\n")
