#!/usr/bin/env python

### Creates plots to compare samples in the treatment groups 
### Takes in .h5ad files from QC_per_sample.py
import matplotlib.pyplot as plt
import scanpy
import pandas as pd
import os
import glob
import warnings
import argparse
warnings.filterwarnings("ignore")
from collections import Counter

# Define the parser
parser = argparse.ArgumentParser(description='Group QC Plots')

parser.add_argument('--in_files', action = "store", dest = "in_files", nargs='*', required=True)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)


# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("Input Files: ", args.in_files)
print("Output folder: ", args.out_dir)
print("***********************")

# Creates a series of plots that compare composition between samples in the groups
def qc_plots_across_groups(adata, out_dir):
    outcome_group_list = adata.obs["condition"].to_list()
    outcome_group_list = list(set(outcome_group_list))

    plt.switch_backend('agg')

    # Creates set of plots for each outcome group (ex: Moderate Covid)
    for outcome in outcome_group_list:
        print("\n"+"Sample Group: ", outcome)

        adata_subset = adata[adata.obs.condition == outcome]

        # Cells per sample
        sample_id_list = adata_subset.obs["sample_id"].to_list()
        counter_obj = Counter(sample_id_list)
        x_keys = list(counter_obj.keys())
        x_values = list(counter_obj.values())
        fig = plt.figure()
        plt.xlabel('Number of Cells Per Sample')
        plt.bar(x_keys, x_values, color ='maroon')
        plt.savefig(os.path.join(out_dir, outcome + "_cells_per_sample.png"))


        # Fraction of Cells per Group
        x_avg = []
        for i in x_values:
            x_avg.append(i/len(sample_id_list))

        fig = plt.figure()
        plt.bar(x_keys, x_avg, color ='purple')
        plt.xlabel('Percent Cells Per Sample')
        plt.savefig(os.path.join(out_dir,outcome + "_percent_cells_per_sample.png"))

        df = pd.DataFrame({'Sample List': x_keys, 'Cells per Sample': x_values, 'Fraction Cells per Sample': x_avg})
        df = df.set_index("Sample List")
        print(df)
        df.to_csv(out_dir+outcome+'_cell_count_summary.csv') 

    # QC plots comparing all samples
    with plt.rc_context():
        scanpy.pl.violin(adata, 'total_counts', jitter=0.4, groupby='sample_id',  multi_panel=True, show=False, rotation=90) 
        plt.savefig(os.path.join(out_dir,"counts_by_sample.png"))

    with plt.rc_context():
        scanpy.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, groupby='sample_id',  multi_panel=True, show=False, rotation=90)
        plt.savefig(os.path.join(out_dir, "number_genes_by_sample.png"))

    with plt.rc_context():
        scanpy.pl.violin(adata, 'pct_counts_mt', jitter=0.4, groupby='sample_id',  multi_panel=True, show=False, rotation=90)
        plt.savefig(os.path.join(out_dir, "percent_mito_by_sample.png"))


# Get a list of files that match the pattern
matched_files = args.in_files

# Loads sample into anndata
adata_list = []
for sample in matched_files:
    print("Loading sample... "+sample)
    adata_list.append(scanpy.read_h5ad(sample))

# Combines sample anndatas into 1 anndata
adata = scanpy.concat(adata_list, join="inner")

print(adata)
qc_plots_across_groups(adata, args.out_dir)

