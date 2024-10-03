#!/usr/bin/env python
# Runs differential expression between 2 conditions including psuedobulk and cell types

import os
import numpy as np
import matplotlib.pyplot as plt
import anndata
import scanpy #as sc
import csv
# Only pandas >= v0.25.0 supports column names with spaces in querys
import warnings
warnings.filterwarnings("ignore")  # plotnine has a lot of MatplotlibDeprecationWarning's
import seaborn as sns
from bioinfokit import visuz
import pandas as pd
from anndata import AnnData
import anndata as ad
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import argparse


# Define the parser
parser = argparse.ArgumentParser(description='pyDeSeq2')

parser.add_argument('--p_adj', dest='p_adj', type=float, default=0.05)
parser.add_argument('--log2fc', dest='log2fc', type=float, default=0.5)
parser.add_argument('--heatmap_num_sig_genes', dest='heatmap_num_sig_genes', type=int, default=50)
parser.add_argument('--min_cells_per_group', dest='min_cells', type=int, default=100)
parser.add_argument('--in_file', dest="in_file", required=True)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("P-value adjusted significance level: ", args.p_adj)
print("Log2 foldchange significance level (absolute value): ", args.log2fc)
print("Number of significant genes for heatmap: ", args.heatmap_num_sig_genes)
print("Minimum number of cells for cell type comparisons: ", args.min_cells)
print("Input file: ", args.in_file)
print("Output folder: ", args.out_dir)
print("***********************")



#Function converts an anndata to be 1 row per sample with the sum of the cells
def sum_by_sample(cell_subset):
    pbs = []
    for sample_id in cell_subset.obs.sample_id.unique():
        sample_cell_subset = cell_subset[cell_subset.obs["sample_id"] == sample_id]
        rep_adata = scanpy.AnnData(X = sample_cell_subset.X.sum(axis = 0), var = sample_cell_subset.var[[]])
        rep_adata.obs_names = [sample_id]
        rep_adata.obs['condition'] = sample_cell_subset.obs['condition'].iloc[0]
        pbs.append(rep_adata)
    cell_subset = scanpy.concat(pbs)
    return cell_subset


def run_pyDeSeq2(adata_group, condition, condition2, cell_subset, cell_subset_2, cell_group, min_cells, out_dir):

    if cell_group == "Psuedobulk":
        cell_subset = sum_by_sample(cell_subset)

        cell_subset_2 = sum_by_sample(cell_subset_2)
    # Further subsets the data by cell type
    else:
        cell_subset = cell_subset[cell_subset.obs.leiden == cell_group]
        print(cell_subset)
        if len(cell_subset.obs_names) < min_cells:
            print("Not enough cells, moving to next type ...")
            return
        cell_subset = sum_by_sample(cell_subset)

        cell_subset_2 = cell_subset_2[cell_subset_2.obs.leiden == cell_group]    
        print(cell_subset_2)
        if len(cell_subset_2.obs_names) < min_cells:
            print("Not enough cells, moving to next type ...")
            return
        cell_subset_2 = sum_by_sample(cell_subset_2)



    #Creates an average counts df for later file write
    cell_sub_df = cell_subset.to_df().transpose()
    cell_sub_df['mean'] = cell_sub_df.mean(axis=1)
    cell_sub_condition_dict = cell_sub_df[['mean']].to_dict()

    cell_sub_df = cell_subset_2.to_df().transpose()
    cell_sub_df['mean'] = cell_sub_df.mean(axis=1)
    cell_sub_condition2_dict = cell_sub_df[['mean']].to_dict()

    cell_subset = ad.concat([cell_subset, cell_subset_2], join="outer")



    # Prep data for pydeseq
    counts = cell_subset.to_df()
    counts = counts.fillna(0) #deseq2 cannot handle na
    counts = counts.round(0).astype(int) #deseq2 only allows intergers
    print(counts)


    # Create deseq object
    dds = DeseqDataSet( 
            counts=counts,
            metadata=cell_subset.obs,
            design_factors="condition")

    # Filtering out genes only found in one sample
    scanpy.pp.filter_genes(dds, min_cells=1) 

    # Runs deseq
    print("\nRunning pyDeSeq2... \n")
    dds.deseq2()

    # Plotting PCA based on condition
    scanpy.tl.pca(dds)
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.pca(dds, color="condition", size=200)
        plt.savefig(os.path.join(out_dir,cell_group+"_"+condition+"_"+condition2+"_Condition_PCA.png"))


    # NOTE: DeseqStats replaces "_" with "-" so you need to change the name of the contrasts if needed
    stat_res = DeseqStats(dds, contrast=["condition", condition, condition2])
    print("PyDeSeq2 Summary\n")
    stat_res.summary()

    # Saves summary info to df
    df = stat_res.results_df

    df.sort_values('padj').to_csv(os.path.join(out_dir,cell_group+"_"+condition+"_"+condition2+"_all_genes_differential_gene_expr_results.csv"))

    # Finds signifcant genes that also had a log 2 fold change > user_input
    sigs = df[(df.padj < args.p_adj) & (abs(df.log2FoldChange) > args.log2fc)]

    print("\n",sigs[sigs.columns[0]].count(), " genes were found to be signficant")

    # Creates file with the averages for each condition
    with open(os.path.join(out_dir,cell_group+"_"+condition+"_"+condition2+"_significant_averages.csv"), "w") as file:
        writer = csv.writer(file)
        writer.writerow([condition+"Average Expr Non-Normalized",condition2+"Average Expr Non-Normalized",'log2FoldChange', 'padj'])
        for index, row in sigs.iterrows():
            writer.writerow([index,cell_sub_condition_dict.get('mean').get(index),cell_sub_condition2_dict.get('mean').get(index),row['log2FoldChange'], row['padj']])
        
    full_sigs = sigs
    # Decreased number of signficant genes for heatmap if there is too many
    if sigs[sigs.columns[0]].count() > args.heatmap_num_sig_genes:
        print("Over ",args.heatmap_num_sig_genes," signifant genes, only using the top ",args.heatmap_num_sig_genes," for the heatmap (by smallest padj)")
        sigs = sigs.sort_values('padj').head(args.heatmap_num_sig_genes)
    
    print(sigs)

    # Creates heatmap if there is more than 1 signifant gene
    if sigs[sigs.columns[0]].count() > 1:
        print("\nCreating the heatmap...")

        # Normalizing and setting up counts for heatmap and volcano plot
        dds.layers['normed_counts']
        dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])

        dds_sigs = dds[:, sigs.index]
        dds_sigs

        grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                            index=dds_sigs.var_names, columns=dds_sigs.obs_names)

        grapher.to_csv(os.path.join(out_dir,cell_group+"_"+condition+"_"+condition2+"_log_dif_expression_table.csv"), sep=',', encoding='utf-8')


        # Creates heatmap
        fig, ax = plt.subplots(figsize=(30, 10))
        sns.clustermap(grapher, z_score=0, cmap = 'RdYlBu_r', yticklabels=True)
        plt.savefig(os.path.join(out_dir,cell_group+"_"+condition+"_"+condition2+"_heatmap_all_top.pdf"), format="pdf", dpi=1000)


    # Creates volcano if there is at least 1 signifant gene down and up regulated
    if (full_sigs["log2FoldChange"] > args.log2fc).any() and (full_sigs["log2FoldChange"] < -args.log2fc).any():


        # Creates volcano plot
        print("\nCreating the volcano plot... ")
        df["padj"] = df["padj"].fillna(1)
        visuz.GeneExpression.volcano(df=df, lfc="log2FoldChange", 
                            pv="padj", 
                            axtickfontname="DejaVu Sans", 
                            axlabelfontname="DejaVu Sans", 
                            lfc_thr=(args.log2fc,args.log2fc),
                            pv_thr=(args.p_adj,args.p_adj),
                            color=("red","grey","blue"), 
                            figname=os.path.join(out_dir,cell_group+"_"+condition+"_"+condition2+"_volcano_plot"), figtype="pdf")

plt.switch_backend('agg')
adata_group = anndata.read_h5ad(args.in_file)
print(adata_group)

#Resets X to be the counts that are not normalized 
adata_group.X = adata_group.layers["counts"]


conditions = list(set(adata_group.obs.condition.to_list()))
conditions.sort()

cell_groups = list(set(adata_group.obs.leiden.to_list()))
cell_groups.append("Psuedobulk")
cell_groups.sort()


count = 1 
print("start of loop")
for condition in conditions:
    for condition2 in conditions[count:]: # Increases start of list to prevent duplicates
        if condition != condition2:
            print("\nProcessing... ", condition, " vs ", condition2, "\n")
            cell_subset = adata_group[adata_group.obs.condition == condition]
            cell_subset_2 = adata_group[adata_group.obs.condition == condition2]
            for cell_group in cell_groups:
                print("For cell type: ", cell_group)
                run_pyDeSeq2(adata_group, condition, condition2, cell_subset, cell_subset_2, cell_group, args.min_cells, args.out_dir)
    count+=1
