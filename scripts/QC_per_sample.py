#!/usr/bin/env python

### Individual sample QC
### Including removal of ambient RNA, doublets, cell/gene/mitochondrial content filtering 
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import scanpy
import gzip
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.interpolate import interpn
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.packages import STAP
import argparse
import warnings
warnings.filterwarnings("ignore")

#R code setup
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

# Define the parser
parser = argparse.ArgumentParser(description='Per sample QC')

parser.add_argument('--sample_name', action="store", dest='sample_name', required=True)
parser.add_argument('--condition', action="store", dest='condition', required=True)
parser.add_argument('--batch', action="store", dest='batch', required=True)
parser.add_argument('--h5_file_name', action="store", dest='h5_file_name', default=0)
parser.add_argument('--h5_raw_file_name', action="store", dest='h5_raw_file_name', default=0)
parser.add_argument('--mtx_tsv_file_name', dest='mtx_tsv_file_name', nargs='*', default=0)
parser.add_argument('--mtx_tsv_raw_file_name', dest='mtx_tsv_raw_file_name', nargs='*', default=0)
parser.add_argument('--cluster_resolution', dest='cluster_resolution', type=float, default=0.4)
parser.add_argument('--min_cells_per_gene', dest='min_cells_per_gene', type=int, default=3)
parser.add_argument('--min_genes_per_cell', dest='min_genes_per_cell', type=int, default=200)
parser.add_argument('--mitochondrial_content_max', dest='mitochondrial_content_max', type=int, default=20)
parser.add_argument('--ribosomal_genelist', dest='ribosomal_genelist', required=True)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("Sample name: ", args.sample_name)
print("Condition: ", args.condition)
print("Batch:", args.batch)
print("Filtered h5 file name:", args.h5_file_name)
print("Raw h5 file name:", args.h5_raw_file_name)
print("Filtered mtx and tsv file list:", args.mtx_tsv_file_name)
print("Raw mtx and tsv file list:", args.mtx_tsv_raw_file_name)
print("Leiden cluster resolution:", args.cluster_resolution)
print("Minimum cells per gene:", args.min_cells_per_gene)
print("Minimum genes per cell:", args.min_genes_per_cell)
print("Mitochondrial content maximum:", args.mitochondrial_content_max)
print("Output folder: ", args.out_dir)
print("***********************")


#converts a matrix to a list
def flatten(matrix):
    return [str(item) for row in matrix for item in row]


# Converts .h5 or mtx/tsv file set into an anndata object
def file_to_adata(sample, condition, batch, file_name, file_type):
    print("Loading ... ", file_name)

    obs_list = [] #holds the sample ids
    batch_list = [] #hold the date the samples were run
    condition_list = [] #hold the treatment condition 

    # Loads data depending on file type
    if file_type.strip() == "h5":
        adata = scanpy.read_10x_h5(file_name)
    elif file_type.strip() == "mtx-tsv":
        for file in file_name:
            if file.strip().endswith("matrix.mtx") or file.strip().endswith("matrix.mtx.gz"):
                adata = scanpy.read_mtx(file).T
            elif file.strip().endswith("barcodes.tsv.gz"):
                barcodes = [line.strip() for line in gzip.open(file, "rt")]
            elif file.endswith("barcodes.tsv"):
                barcodes = [line.strip() for line in open(file, "r")]
            elif file.endswith("genes.tsv.gz"):
                genes = [line.strip() for line in gzip.open(file, "rt")]
            elif file.endswith("genes.tsv"):
                genes = [line.strip() for line in open(file, "r")]

        adata.obs_names = barcodes
        adata.var_names = genes


    #Adds additional annotation to obs
    obs_list.append([sample] * len(adata.obs))
    condition_list.append([condition] * len(adata.obs))
    batch_list.append([batch] * len(adata.obs))
    obs_list = flatten(obs_list)
    condition_list = flatten(condition_list)
    batch_list = flatten(batch_list)
    adata.obs["sample_id"] = obs_list
    adata.obs["condition"] = condition_list
    adata.obs["batch"] = batch_list

    adata.var_names_make_unique()

    return(adata)



# Basic QC (removal of ambient RNA, doublets, cell/gene/mitochondrial content filtering)
# Saves results in .h5ad file
def quaility_analysis(adata, sample, raw_file_name, raw_file_name_type, cluster_res, min_cells_per_gene, min_genes_per_cell, mito_content_max, ribosomal_genelist, out_dir):
    prefix = sample
    print("\n Orginal Filtered Ouput from Cell Ranger")
    print(adata)
    plt.switch_backend('agg')

    # Perform SVD 
    tsvd = TruncatedSVD(n_components=2)
    tsvd.fit(adata.X)
    X = tsvd.transform(adata.X)

    # Plot the cells in the 2D PCA projection
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(X[:,0], X[:,1], alpha=0.5, c="green")
    plt.savefig(os.path.join(out_dir, prefix + "_2D_PCA_projection.png"))


    # Test for library saturation
    # Create a plot showing genes detected as a function of UMI counts.
    fig, ax = plt.subplots(figsize=(10, 7))
    x = np.asarray(adata.X.sum(axis=1))[:,0]
    y = np.asarray(np.sum(adata.X>0, axis=1))[:,0]
    ax.scatter(x, y, color="green", alpha=0.25)
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.savefig(os.path.join(out_dir, prefix + "_library_saturation.png"))


    # Creates knee plot
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    cell_set = np.arange(len(knee))
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(knee, cell_set, linewidth=5, color="g")
    ax.set_xlabel("UMI Counts")
    ax.set_ylabel("Set of Barcodes")
    plt.grid(True, which="both")
    plt.savefig(os.path.join(out_dir, prefix + "_knee_plot.png"))
    

    # Checks if raw samples were provided so ambient RNA can be removed
    if raw_file_name != 0: 
        print("\nRemoving ambient RNA ...."+"\n")
        adata = cook_soup(adata, raw_file_name, raw_file_name_type)
        print(adata)

    # Creates a copy of the anndata obj (the orignal anndata gets normalized for doublet plotting)
    whole_adata = adata.copy()
    # whole_adata.var_names_make_unique()


    #Uses Scrublet to remove doublets
    print("\nRemoving doublets ...")
    scanpy.external.pp.scrublet(adata)
    print(adata.obs["predicted_doublet"].value_counts())
    adata.obs["predicted_doublet_category"] = adata.obs["predicted_doublet"].astype("category")
    doublets_df = adata.obs["predicted_doublet"].to_frame()
    doublets_df = doublets_df[doublets_df.predicted_doublet == True]

    #Visuallize doublet results
    scanpy.pp.normalize_total(adata, target_sum = 1e4)
    scanpy.pp.log1p(adata)
    scanpy.tl.pca(adata)
    scanpy.pp.neighbors(adata)
    scanpy.tl.umap(adata)
    scanpy.tl.leiden(adata, resolution = cluster_res)

    with plt.rc_context({'figure.figsize': (4, 4)}):
        scanpy.external.pl.scrublet_score_distribution(adata)
        plt.savefig(os.path.join(out_dir, prefix + "_scrublet_score_distribution.png"))

    with plt.rc_context({'figure.figsize': (4, 4)}):
        scanpy.pl.umap(adata, color = ['leiden', 'predicted_doublet_category'])
        plt.savefig(os.path.join(out_dir, prefix + "_leiden_doublet.png"))


    # Remove doublets from non normalized anndata
    whole_adata.obs['doublet'] = whole_adata.obs.index.isin(doublets_df.index)
    adata = whole_adata[~whole_adata.obs.doublet]


    # Calculates mitochondrial and ribosomal content
    print("\nPerforming count based filtering ...")
    adata.var['mt'] = adata.var.index.str.startswith('MT-')
    #ribosomal gene list from: http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt
    #ribo_url = "reference_files/KEGG_RIBOSOME.v2024.1.Hs.txt"
    ribo_genes = pd.read_table(ribosomal_genelist, skiprows=2, header = None)
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    scanpy.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

    # Filters data based on gene and cell counts
    scanpy.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    scanpy.pp.filter_cells(adata, min_genes=min_genes_per_cell)

    # Creates count distribution plots
    adata.var.sort_values('n_cells_by_counts')
    adata.obs.sort_values('n_genes_by_counts')
    with plt.rc_context():
        scanpy.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], jitter=0.4, multi_panel=True)
        plt.savefig(os.path.join(out_dir, prefix + "_count_distributions.png"))


    # Removes mitochondrial content over a certain number (ex: 20, keeps every cell < 20%)
    adata = adata[adata.obs.pct_counts_mt < mito_content_max]
    with plt.rc_context():
        scanpy.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
        plt.savefig(os.path.join(out_dir, prefix + "_percent_mito_by_total_counts.png"))

    # Save filtered data 
    print("\nSaving filtered data...")
    print("Located at: ", os.path.join(out_dir, prefix + "_qc_filtered.h5ad") )
    adata_write = adata.copy()
    adata_write.write_h5ad(os.path.join(out_dir, prefix + "_qc_filtered.h5ad"))

    return adata
    



# Creates groups for soupX
# SoupX code is adpated from mousepixels: soupX/soupX_python_test.ipynb
def get_soupx_group(adata):
    adata_pp = adata.copy()
    scanpy.pp.normalize_per_cell(adata_pp)
    scanpy.pp.log1p(adata_pp)
    scanpy.pp.pca(adata_pp)
    scanpy.pp.neighbors(adata_pp)
    scanpy.tl.leiden(adata_pp, key_added="soupx_groups")
    return adata_pp.obs['soupx_groups']
    
    
# Prepare parameters for SoupX
def prepare_broth(adata, raw_file_name, raw_file_name_type):

    #get raw data, puts placeholders for the batch, condition, sample name as they are not used in downstream analysis
    raw_adata = file_to_adata("raw",0,0, raw_file_name, raw_file_name_type)

    #Removes genes only found in raw or filtered
    common_genes = adata.var_names.intersection(raw_adata.var_names)
    adata = adata[:, common_genes]
    raw_adata = raw_adata[:, common_genes]

    # Make into individual components to pass to R
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    raw_adata = raw_adata.X.T
    
    #get leiden clusters
    soupx_groups = get_soupx_group(adata)

    return data, adata, raw_adata, genes, cells, soupx_groups

# Defines SoupX R code
r_soupX_func = """    
library(SoupX)

make_soup <- function(data, raw, genes, cells, soupx_groups){
    # specify row and column names of data
    rownames(data) = genes
    colnames(data) = cells
    # ensure correct sparse format for table of counts and table of droplets
    data <- as(data, "sparseMatrix")
    raw <- as(raw, "sparseMatrix")

    # Generate SoupChannel Object for SoupX 
    sc = SoupChannel(raw, data, calcSoupProfile = FALSE)

    # Add extra meta data to the SoupChannel object
    soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
    sc = setSoupProfile(sc, soupProf)
    # Set cluster information in SoupChannel
    sc = setClusters(sc, soupx_groups)

    # Estimate contamination fraction
    sc  = autoEstCont(sc, doPlot=FALSE)
    # Infer corrected table of counts and round to integer
    out = adjustCounts(sc, roundToInt = TRUE)
    
    return(out)
}
"""
r_soupX = STAP(r_soupX_func, "r_pkg")

# Calls SoupX functions
def cook_soup(adata, raw_file_name, raw_file_name_type):

    #Sets up SoupX paramaters
    data, adata, raw, genes, cells, soupx_groups = prepare_broth(adata, raw_file_name, raw_file_name_type)

    # Execute the R code and get the corrected counts
    out = r_soupX.make_soup(data, raw, genes, cells, soupx_groups)

    adata.layers["raw_counts"] = adata.X
    adata.layers["soupX_counts"] = out.T
    adata.X = adata.layers["soupX_counts"]
    
    return adata

# Processing sample driver function
def process_sample(sample, condition, batch, file_name, file_name_type, raw_file_name, raw_file_name_type, cluster_res, min_cells_per_gene, min_genes_per_cell, mito_content_max, ribosomal_genelist, out_dir):
    print("\n Processing sample ", sample)

    adata = file_to_adata(sample, condition, batch, file_name, file_name_type)

    adata = quaility_analysis(adata, sample, raw_file_name, raw_file_name_type, cluster_res, min_cells_per_gene, min_genes_per_cell, mito_content_max, ribosomal_genelist, out_dir)
    
    print("Final anndata")
    print(adata)



# Checks that there is a valid configuration of files and identifies file types
def set_files(args):
    if args.h5_file_name != 0:
        file_name = args.h5_file_name
        file_type = "h5"
    elif isinstance(args.mtx_tsv_file_name, list) and len(args.mtx_tsv_file_name) == 3:
        file_name = args.mtx_tsv_file_name
        file_type = "mtx-tsv"   
    else:
        raise ValueError("No file present or not a valid type. Please use the cell ranger output filtered h5 file or the filtered mtx-tsv files.")
    
    if args.h5_raw_file_name != 0:
        raw_file_name = args.h5_raw_file_name
        raw_file_type = "h5"
    elif  isinstance(args.mtx_tsv_raw_file_name, list) and len(args.mtx_tsv_raw_file_name) == 3:
        raw_file_name = args.mtx_tsv_raw_file_name
        raw_file_type = "mtx-tsv" 
    else:
        raw_file_name = 0
        raw_file_type = 0       
    
    return file_name, file_type, raw_file_name, raw_file_type

try:
    file_name, file_type, raw_file_name, raw_file_type = set_files(args)
except Exception as e:
    print(e)
    sys.exit()


process_sample(args.sample_name, args.condition, args.batch, file_name, file_type, 
               raw_file_name, raw_file_type, args.cluster_resolution, args.min_cells_per_gene, 
               args.min_genes_per_cell, args.mitochondrial_content_max, args.ribosomal_genelist, args.out_dir)
