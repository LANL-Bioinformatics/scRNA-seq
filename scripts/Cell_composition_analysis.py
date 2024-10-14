#!/usr/bin/env python
# Plots percent of each cell type and runs scCODA to deterine credible differences between condition groups

import matplotlib.pyplot as plt
import anndata
# Only pandas >= v0.25.0 supports column names with spaces in querys
import warnings
warnings.filterwarnings("ignore")  # plotnine has a lot of MatplotlibDeprecationWarning's
import pandas as pd
import pertpy as pt #contains the scCODA library
import argparse
import os

# Define the parser
parser = argparse.ArgumentParser(description='Cell Composition Analysis')

parser.add_argument('--p_adj', dest='p_adj', type=float, default=0.05)
parser.add_argument('--in_file', dest='in_file', required=True)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("P-value adjusted significance level: ", args.p_adj)
print("Input file: ", args.in_file)
print("Output folder: ", args.out_dir)
print("***********************")


# Read data file with cell type names
out_dir = args.out_dir
adata = anndata.read(args.in_file)

# Uses colors from leiden clusters for barplots
color_set = adata.uns["leiden_colors"]


df = pd.DataFrame(columns=adata.obs["merged leiden"].unique())

# Creates df with the columns as the different cell types and the rows are the outcomes with the data values being the percent cells (decimal)
for con in adata.obs.condition.unique():
    adata_subset = adata[adata.obs["condition"] == con]
    total_cells = len(adata_subset.obs_names)
    groups_avg = []
    for group_num in adata.obs["merged leiden"].unique():
        adata_subset_group = adata_subset[adata_subset.obs["merged leiden"] == group_num]
        groups_avg.append(len(adata_subset_group.obs_names)/total_cells)
    df.loc[con] = groups_avg

df.to_csv(out_dir+"composition_percent_by_cell_type.csv")
print(df)

# Create barplot of cell type percentages seperated out by outcome
plt.switch_backend('agg')
with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    df.plot.barh(stacked=True, color=color_set)
    plt.legend(loc='center right',bbox_to_anchor=(1.5,0.5)) 
    plt.savefig(os.path.join(out_dir,"cell_type_by_condition_composition_bar_plot.png"), bbox_inches='tight')



# resets X to non normalized values for scCODA
adata.X = adata.layers["counts"]

# Creates a boxplot with all the outcomes
sccoda = pt.tl.Sccoda()
mdata = sccoda.load(adata, type="cell_level", generate_sample_level=True, cell_type_identifier="leiden", sample_identifier="sample_id", covariate_obs=["condition"])
with plt.rc_context():  
    fig, ax = plt.subplots(figsize=(12, 7))
    sccoda.plot_boxplots(mdata, feature_name="condition")
    plt.savefig(os.path.join(out_dir,"whole_dataset_coda_boxplot.png"), bbox_inches='tight')


# Function runs scCODA between 2 conditions in adata
def run_sccoda(condition):

    # Creates instance and aggreates data
    sccoda_model = pt.tools.Sccoda()
    sccoda_data = sccoda_model.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier="leiden",
        sample_identifier="sample_id",
        covariate_obs=["condition"],
    )

    # Select the 2 data types (outcomes)
    sccoda_data.mod["coda_mod"] = sccoda_data["coda"][
        sccoda_data["coda"].obs["condition"].isin([condition[0], condition[1]])
    ].copy()
    print(sccoda_data["coda_mod"])

    # Creates a boxplot of the 2 conditions
    with plt.rc_context():
        sccoda.plot_boxplots(sccoda_data, modality_key="coda_mod", feature_name="condition", add_dots=True)
        plt.savefig(os.path.join(out_dir, condition[0]+"_"+condition[1]+"_coda_boxplot.png"), bbox_inches='tight')

    # Model setup and inference
    sccoda_data = sccoda_model.prepare(
        sccoda_data,
        modality_key="coda_mod",
        formula="condition",
        reference_cell_type="automatic",
    )
    sccoda_data["coda_mod"]

    # Run MCMC
    sccoda_model.run_nuts(sccoda_data, modality_key="coda_mod",rng_key=1234)
    sccoda_data["coda_mod"]

    # Sets FDR to 0.05
    print(sccoda_model.set_fdr(sccoda_data, modality_key="coda_mod", est_fdr=args.p_adj))

    print(sccoda_model.summary(sccoda_data, modality_key="coda_mod"))

    print(sccoda_model.credible_effects(sccoda_data, modality_key="coda_mod"))

    # Creates a pandas series with the status (True/False)
    credible = sccoda_model.credible_effects(sccoda_data, modality_key="coda_mod")
  
    # Checks if there were credible differences found
    if True in credible.values:
        print("At least 1 difference was found to be credible, creating a plot....")

        # creates a barplot of the credible results
        with plt.rc_context():
            sccoda.plot_effects_barplot(sccoda_data, modality_key="coda_mod", parameter="Final Parameter")
            plt.savefig(os.path.join(out_dir, condition[0]+"_"+condition[1]+"_credible_coda_bar_plot.png"), bbox_inches='tight')


conditions = list(set(adata.obs.condition.to_list()))
conditions.sort()

count = 1 
for condition in conditions:
    for condition2 in conditions[count:]: # Increases start of list to prevent duplicates
        if condition != condition2:
            print("\nProcessing... ", condition, " vs ", condition2, "\n")
            run_sccoda([condition, condition2])

    count+=1

