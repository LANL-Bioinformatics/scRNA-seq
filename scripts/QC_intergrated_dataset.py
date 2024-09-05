#adpated from juypter notebook
#added additional analysis from sambomics video and put all treatment groups together 
#Removed clustering and cell typing so it is just QC
#Changed

import matplotlib.pyplot as plt
import scanpy #as sc
import os
import glob
from scipy.sparse import csr_matrix
# from scipy.io import mmread
import pandas as pd
import argparse
# Only pandas >= v0.25.0 supports column names with spaces in querys
import warnings
warnings.filterwarnings("ignore") 



# Define the parser
parser = argparse.ArgumentParser(description='Intergrated Dataset QC')

parser.add_argument('--HVG_num_top_genes', dest='HVG_num_top_genes', type=int, default=5000)
parser.add_argument('--PCA_N_Comps', dest='PCA_N_Comps', type=int, default=20)
parser.add_argument('--Neighbors_Number', dest='Neighbors_Number', type=int, default=30)
parser.add_argument('--Neighbors_N_PCS', dest='Neighbors_N_PCS', type=int, default=20)
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("Highly Varible Genes number top genes: ", args.HVG_num_top_genes)
print("PCA number of components: ", args.PCA_N_Comps)
print("Number of neighbors: ", args.Neighbors_Number)
print("Neighbors number PCS: ", args.Neighbors_N_PCS)
print("Output folder: ", args.out_dir)
print("***********************")


# #converts a matrix to a list
# def flatten(matrix):
#     return [str(item) for row in matrix for item in row]

#Turns tuples into a dictionary
def Convert(tup, di):
    for a, b in tup:
        di[a] = b
    return di

def intergrated_dataset_qc(adata, HVG_num_top_genes, PCA_N_Comps, Neighbors_Number, Neighbors_N_PCS, out_dir):
    print(adata)
    adata.X = csr_matrix(adata.X)

    adata.layers["counts"] = adata.X.copy()

    #Normalization
    print("\nNormalizing data... ")
    scanpy.pp.normalize_total(adata, target_sum=1e4)
    scanpy.pp.log1p(adata)

    adata.raw = adata

    #Umap plot for samples
    print("\nApplying batch correction... ")
    scanpy.pp.highly_variable_genes(adata, n_top_genes=HVG_num_top_genes)
    scanpy.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, n_comps=PCA_N_Comps)
    scanpy.pp.neighbors(adata, n_neighbors=Neighbors_Number, n_pcs=Neighbors_N_PCS)
    scanpy.tl.umap(adata)
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='sample_id', ax=ax)
        plt.savefig(out_dir+"sample_id_clustering.png", bbox_inches='tight')

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='batch', ax=ax)
        plt.savefig(out_dir+"batch_clustering.png", bbox_inches='tight')

    scanpy.external.pp.harmony_integrate(adata, 'sample_id', adjusted_basis='X_pca', max_iter_harmony=50)
    print(adata)
    scanpy.pp.neighbors(adata, n_neighbors=Neighbors_Number, n_pcs=Neighbors_N_PCS)
    scanpy.tl.umap(adata)

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='sample_id', ax=ax)
        plt.savefig(out_dir+"sample_id_clustering_batch_effect_corrected.png",bbox_inches='tight')

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.umap(adata, color='batch', ax=ax)
        plt.savefig(out_dir+"batch_clustering_batch_effect_corrected.png",bbox_inches='tight')


    adata.obs['condition'] = adata.obs['condition'].astype('category')
    sample_means = pd.DataFrame(adata.to_df()).groupby(adata.obs['condition']).mean()
    adata_sample = scanpy.AnnData(X=sample_means.values)
    adata_sample.obs['condition'] = sample_means.index
    print(adata_sample)
    scanpy.tl.pca(adata_sample)

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.pca(adata_sample, color='condition', ax=ax, size=1)
        plt.savefig(out_dir+"PCA_1_PCA_2.png",bbox_inches='tight')

    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(10, 7))
        scanpy.pl.pca(adata, color="condition", size=200)
        plt.savefig(out_dir+"Condition_PCA.png")


    adata.write_h5ad(out_dir+'integrated.h5ad')


# Get a list of files that match the pattern
matched_files = glob.glob(os.path.join(args.out_dir, "*_qc_filtered.h5ad"))

# Loads sample into anndata
adata_list = []
for sample in matched_files:
    print("Loading sample... "+sample)
    adata_list.append(scanpy.read_h5ad(sample))


# Combines sample anndata into 1 anndata
adata = scanpy.concat(adata_list, join="outer")

intergrated_dataset_qc(adata, args.HVG_num_top_genes, args.PCA_N_Comps, args.Neighbors_Number, args.Neighbors_N_PCS, args.out_dir)