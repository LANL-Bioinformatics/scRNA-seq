#!/usr/bin/env python

# Normalizes data, runs Harmony batch correction and saves dataset in integrated.h5ad
# Takes in input .h5ad file from QC_per_sample.py

import matplotlib.pyplot as plt
import scanpy
import os
import glob
from scipy.sparse import csr_matrix
import argparse
import warnings
warnings.filterwarnings("ignore") 
import pandas as pd


# Define the parser
parser = argparse.ArgumentParser(description='Intergrated Dataset QC')

parser.add_argument('--HVG_num_top_genes', dest='HVG_num_top_genes', type=int, default=5000)
parser.add_argument('--PCA_N_Comps', dest='PCA_N_Comps', type=int, default=30)
parser.add_argument('--Neighbors_Number', dest='Neighbors_Number', type=int, default=30)
parser.add_argument('--Neighbors_N_PCS', dest='Neighbors_N_PCS', type=int, default=20)
parser.add_argument('--Scale_Max', dest='Scale_Max', type=int, default=20)
parser.add_argument('--in_files', dest='in_files', nargs="*", required=True)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("Highly Varible Genes number top genes: ", args.HVG_num_top_genes)
print("PCA number of components: ", args.PCA_N_Comps)
print("Number of neighbors: ", args.Neighbors_Number)
print("Neighbors number PCS: ", args.Neighbors_N_PCS)
print("Scale data to max value:  ", args.Scale_Max)
print("Input files: ", args.in_files)
print("Output folder: ", args.out_dir)
print("***********************")


# Normalizes and batch corrects
def intergrated_dataset_qc(adata, HVG_num_top_genes, PCA_N_Comps, Neighbors_Number, Neighbors_N_PCS, Scale_Max, out_dir):
    print(adata)

    adata.layers["counts"] = adata.X.copy()

    # Normalization
    print("\nNormalizing data... ")
    scanpy.pp.normalize_total(adata, target_sum=1e6)
    scanpy.pp.log1p(adata)
    scanpy.pp.scale(adata, zero_center=True, max_value=Scale_Max)
 

    adata.raw = adata

    # Set up code for batch correction
    print("\nApplying batch correction... ")
    scanpy.pp.highly_variable_genes(adata, n_top_genes=HVG_num_top_genes)
    scanpy.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, n_comps=PCA_N_Comps)
    scanpy.pp.neighbors(adata, n_neighbors=Neighbors_Number, n_pcs=Neighbors_N_PCS)
    scanpy.tl.umap(adata, min_dist=0.3)

    # Plots UMAP without batch correction
    plt.switch_backend('agg')
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='sample_id', ax=ax)
        plt.savefig(os.path.join(out_dir, "sample_id_clustering.png"), bbox_inches='tight')
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='batch', ax=ax)
        plt.savefig(os.path.join(out_dir, "batch_clustering.png"), bbox_inches='tight')

    # Runs Harmony batch correction
    scanpy.external.pp.harmony_integrate(adata, 'sample_id', adjusted_basis='X_pca', max_iter_harmony=50)
    scanpy.pp.neighbors(adata, n_neighbors=Neighbors_Number, n_pcs=Neighbors_N_PCS)
    scanpy.tl.umap(adata, min_dist=0.3)

    # Plots UMAP without batch correction
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='sample_id', ax=ax)
        plt.savefig(os.path.join(out_dir, "sample_id_clustering_batch_effect_corrected.png"),bbox_inches='tight')
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='batch', ax=ax)
        plt.savefig(os.path.join(out_dir, "batch_clustering_batch_effect_corrected.png"),bbox_inches='tight')

    # Plots varaince ratio of PCA components
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.pca_variance_ratio(adata, log=True)
        plt.savefig(os.path.join(out_dir, "variance_ratio_PCA.png"),bbox_inches='tight')

    print("\nWriting dataset to integrated.h5ad ...")
    adata.write_h5ad(os.path.join(out_dir, 'integrated.h5ad'))


# Get a list of files that match the pattern
matched_files = args.in_files

#Loads sample into anndata
adata_list = []
for sample in matched_files:
    print("Loading sample... "+sample)
    adata_list.append(scanpy.read_h5ad(args.out_dir+sample))

# Combines sample anndata into 1 anndata
adata = scanpy.concat(adata_list, join="inner")


intergrated_dataset_qc(adata, args.HVG_num_top_genes, args.PCA_N_Comps, args.Neighbors_Number, args.Neighbors_N_PCS, args.Scale_Max, args.out_dir)
