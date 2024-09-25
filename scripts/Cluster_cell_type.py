#!/usr/bin/env python
# Uses a list of user provided marker genes (or default) to assign cell types to the leiden clusters

import os
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import anndata
import scanpy
import warnings
warnings.filterwarnings("ignore")  # plotnine has a lot of MatplotlibDeprecationWarning's
import pandas as pd
import argparse
import numpy as np


# Define the parser
parser = argparse.ArgumentParser(description='Cluster Cell Type')

parser.add_argument('--cluster_resolution', dest='cluster_resolution', type=float, default=0.4)
parser.add_argument('--percentile', dest='percentile', type=int, default=80)
parser.add_argument('--cell_types', dest='cell_types', nargs='*', default= ["T cells: D3G,TRBC2,CD3D,CD3E,IL7R,LTB", "NK cells:TYROBP,FCGR3A,NKG7,TRDC,KLRF1,KLRD1,GNLY", "B cells: MS4A1,PXK,CD19,CD74,CD79A,IGHM", "Plasma cells: JCHAIN,MZB1,IGHG1,SPAG4", "Proliferating lymphocytes: MKI67,CD3G,FCGR3A", "Monocytes: CD14, FCGR3A,LYZ,CFP,APOBEC3A,CCR2", "cDCs: CD1C,BDCA4,FCER1A", "pDCs: SERPINF1,BST2,MAP3K2,KLK1,TRADD,CLEC4C", "platelets: PF4,GP1BA,SELP,PPBP,ITGA2B", "Erythrocytes: HBA2,HBA1,HBB"])
parser.add_argument('--integrated_data', dest= 'integrated_data', required = True)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("Leiden cluster resolution:", args.cluster_resolution)
print("Gene expression values above this percentile are taken for cell typing: ", args.percentile)
print("Cell Types: ", args.cell_types)
print("Input File: ", args.integrated_data)
print("Output folder: ", args.out_dir)
print("***********************")


# Load data from intergrated sample QC output
out_dir = args.out_dir
adata = anndata.read(args.integrated_data)

# Leiden Clustering 
plt.switch_backend('agg')

scanpy.tl.leiden(adata, resolution = args.cluster_resolution)

with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    scanpy.pl.umap(adata, color='leiden', ax=ax)
    plt.savefig(os.path.join(out_dir, "UMAP_leiden_All.png"))


# Count cells in each Leiden cluster
cluster_counts = adata.obs["leiden"].value_counts()
cluster_counts = cluster_counts.sort_values(ascending=False)

# Calculated scores for each marker gene
df_cell_type_summary = pd.DataFrame()
for item in args.cell_types:
    gene_list = []
    cell_type_group = item.split(":")[0]
    print(cell_type_group)
    cell_group_summary = pd.DataFrame()

    # Iterates through all marker genes assigned to a cell type
    for gene in item.split(":")[1].split(","):
        gene = gene.strip()
        if gene in adata.var_names:

            # Create Data Frame with leiden, cell id and gene column from adata
            df = adata.obs[["leiden"]].copy()
            df["cell_id"] = adata.obs.index
            df[gene] = adata[:, gene].X.flatten().toarray()
            df = df.reset_index(drop=True)

            # Figures out the expression level at the top percentile and creates a dataframe only containing values above it
            print(args.percentile)
            threshold = np.percentile(df[gene], args.percentile)
            highest_expr_data = df.loc[df[gene] > threshold]
            cell_counts = highest_expr_data['leiden'].value_counts()
            
            # Created new df with leiden cluster cell counts, totals cells per cluster and an average
            summary = pd.concat([cell_counts, cluster_counts], axis=1)
            summary.columns = ["cells", "total_cells"]
            summary["Average "+gene] = summary["cells"]/summary["total_cells"]
            
            # Adds average for that marker gene to a new column in the summary df for the cell type
            cell_group_summary = pd.concat([cell_group_summary , summary["Average "+gene]], axis=1)


    # Adds average group expression for the cell type to the cell type summary dataframe and saves marker gene summary df to file
    cell_group_summary[cell_type_group] = cell_group_summary.mean(axis=1)
    print(cell_group_summary)
    print("***********************")

    cell_group_summary.to_csv(os.path.join(out_dir,cell_type_group+"_marker_gene_average.csv"))
    df_cell_type_summary = pd.concat([df_cell_type_summary , cell_group_summary[cell_type_group]], axis=1)


# Identifies the cell type with the highest score for each leiden group
print("\nFull DataSet Summary")
df_cell_type_summary["Top Cell Type"] = df_cell_type_summary.idxmax(axis=1)
df_cell_type_summary.index = pd.to_numeric(df_cell_type_summary.index)
df_cell_type_summary = df_cell_type_summary.sort_index()
print(df_cell_type_summary)
df_cell_type_summary.to_csv(os.path.join(out_dir,"all_types_marker_gene_average.csv"))

calculated_cell_types_gene = df_cell_type_summary["Top Cell Type"].to_numpy()


# Rename groups with calculated names
unique_calculated_cell_types_gene = []
count = 0
for cell_type in calculated_cell_types_gene:
    if cell_type in unique_calculated_cell_types_gene:
        unique_calculated_cell_types_gene.append(cell_type+" "+str(count))
    else:
        unique_calculated_cell_types_gene.append(cell_type)


    count+=1


calculated_cell_types_gene = unique_calculated_cell_types_gene


#Renames groups
print("\nFinal Cell Types ...")
print(calculated_cell_types_gene)
adata.rename_categories('leiden', calculated_cell_types_gene)

#Saves new h5ad file with updated names 

adata.write_h5ad(Path(os.path.join(out_dir, "adata_w_leiden_groups.h5ad")))

#Replots the leiden clustering with new groups
with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    scanpy.pl.umap(adata, color='leiden', ax=ax)
    plt.savefig(os.path.join(out_dir,"UMAP_leiden_automated_names.png"), bbox_inches='tight')
with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    scanpy.pl.umap(adata, color='leiden', legend_loc='on data')
    plt.savefig(os.path.join(out_dir,"UMAP_leiden_automated_names_labeled_on_plot_.png"), bbox_inches='tight')


# Custom colormap
colors = [(0.85, 0.85, 0.85), (1, 0, 0)]  # light grey to red
cmap = LinearSegmentedColormap.from_list("custom_grey_red", colors)


# Creates the same leiden clustering plots as earlier by colored by expression of marker gene list from the user
print("\nCreating gene colored leiden clustering plots ...")
for item in args.cell_types:
    gene_list = []
    cell_type_group = item.split(":")[0]
    for gene in item.split(":")[1].split(","):
        gene = gene.strip()
        # Checks if gene is present in the dataset
        if gene in adata.var_names:
            gene_list.append(gene)

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color=gene_list, color_map=cmap)
        plt.savefig(os.path.join(out_dir,"UMAP_leiden_"+cell_type_group+"_marker_gene.png"), bbox_inches='tight')



# Created UMAP plots for each condition
print("\nCreating UMAP plots for each treatment condition...")
condition_list = list(set(adata.obs["condition"].to_list()))
for condition in condition_list:
    subset_adata = adata[adata.obs.condition == condition]

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(subset_adata, color='leiden', ax=ax)
        plt.savefig(os.path.join(out_dir,"UMAP_leiden_"+condition+"_automated_names.png"), bbox_inches='tight')
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(subset_adata, color='leiden', legend_loc='on data')
        plt.savefig(os.path.join(out_dir,"UMAP_leiden_"+condition+"_automated_names_labeled_on_plot_.png"), bbox_inches='tight')


