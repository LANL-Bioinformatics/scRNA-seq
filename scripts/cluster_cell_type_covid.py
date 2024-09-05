# adpated from juypter notebook cluster_analysis.ipynb
# Reads input from run_qc_qa.py for gene expression data
# Reads input from convert_cell_ranger_h5ad.py for cell surface protein data
# Assigns cell types to leiden groups using diferentially expressed marker genes
# conda env: 10x_py_scvi



# import sys
import numpy as np
import matplotlib.pyplot as plt
import anndata
import scanpy
# Only pandas >= v0.25.0 supports column names with spaces in querys
import warnings
warnings.filterwarnings("ignore")  # plotnine has a lot of MatplotlibDeprecationWarning's
import pandas as pd
import anndata as ad
from collections import Counter
import os
import scvi


### Variable Definitions ###
file_list = "/panfs/biopan04/unr_thailand/unr_analysis/file_lists/GEX_file_list.txt"
out_dir = "/panfs/biopan04/unr_thailand/unr_analysis/figs/analysis_v2/covid/" 

only_gene_expr = True
outcome_group = "All" #Set to All or specific group ie Dead 


### Helper Functions ###
def flatten(matrix):
    return [str(item) for row in matrix for item in row]

def find_maxes(my_dict):
        scores = {}
        high_score = 0
        for key, value in my_dict.items():
            scores.setdefault(value, []).append(key)
            if value > high_score:
                high_score = value
        results = scores[high_score]
        if len(results) == 1:
            return [results[0], high_score]
        else:
            my_str = ""
            for x in results:
                my_str = my_str + "|" + x
            my_str = my_str[1:]
            return [my_str, high_score]
        
def find_score(most_common_cell_type, i):  
    score = 0
    print(most_common_cell_type, " ", calculated_cell_types_gene[i])
    print("Now testing ", most_common_cell_type, " with current score: ", score)
    if calculated_cell_types_gene[i] == most_common_cell_type:
        print("matches at gene with score ", calculated_cell_types_gene_score[i])
        score+=calculated_cell_types_gene_score[i]
    if calculated_cell_types_prot_high[i] == most_common_cell_type:
        score+=calculated_cell_types_prot_high_score[i]
        print("matches at prot high with score ", calculated_cell_types_prot_high_score[i])
    if calculated_cell_types_prot_low[i] == most_common_cell_type:
        score+=calculated_cell_types_prot_low_score[i]
        print("matches at prot low with score ", calculated_cell_types_prot_low_score[i])
    print("The score is now ", score)
    return round(score, 2)

# prot_files = [f for f in os.listdir("/panfs/biopan04/unr_thailand/cell_ranger_CSP_h5ad_new_antibody_list")]
# prot_files.sort()
# print("FILE NAMES     ", prot_files)



### Cluster Analysis Setup ###


### Sorts the protien dataset into treatment groups
with open("/panfs/biopan04/unr_thailand/unr_analysis/file_lists/cell_type_markers_expanded_list.txt") as f:
    cell_type_infos = f.read().splitlines()
    cell_type_infos = cell_type_infos[1:len(cell_type_infos)] #removes header 

cell_types_prot_pos = {}
cell_types_prot_neg = {}
cell_types_gene = {}

for cell_type_info in cell_type_infos:
    item_list = cell_type_info.split("\t")

    #Protien high Expression
    for item in item_list[4].split(","):
        item = item.strip()
        item = item.strip('"')
        if cell_types_prot_pos.get(item) == None:
            cell_types_prot_pos[item] = [item_list[0].strip().strip(',')]
        else:
             cell_types_prot_pos[item].append(item_list[0].strip().strip(','))

    #Protien low Expression
    for item in item_list[5].split(","):
        item = item.strip()
        item = item.strip('"')
        if item == "":
            print("no down regualted")
        elif cell_types_prot_neg.get(item) == None:
            cell_types_prot_neg[item] = [item_list[0].strip().strip(',')]
        else:
             cell_types_prot_neg[item].append(item_list[0].strip().strip(','))


    #RNA Gene Data         
    for item in item_list[3].split(","):
        item = item.strip()
        item = item.strip('"')
        if cell_types_gene.get(item) == None:
            cell_types_gene[item] = [item_list[0].strip().strip(',')]
        else:
             cell_types_gene[item].append(item_list[0].strip().strip(','))




### Load Data Gene ###
adata = anndata.read(out_dir+"integrated.h5ad")


# ### Load Data Cell Surface Protein ###
# prot_files = [f for f in os.listdir("/panfs/biopan04/unr_thailand/cell_ranger_CSP_h5ad_new_antibody_list")]
# prot_files.sort()
# print("FILE NAMES     ", prot_files)


# sample_id_list = []
# cell_id_list = []
# group = prot_files
# print(group)
# first_sample=group[0]

# adata_group_prot = anndata.read("/panfs/biopan04/unr_thailand/cell_ranger_CSP_h5ad_new_antibody_list/"+first_sample)

# print(adata_group_prot.var.gene_id)
# gene_id = adata_group_prot.var.gene_id
# gene_id = gene_id.to_dict()
# gene_id = gene_id.values()
# index_gene_id_key = dict(zip(adata_group_prot.var_names, gene_id))

# print(first_sample.split(".")[0])

# sample_id_list.append([first_sample.split(".")[0]] * len(adata_group_prot.obs))
# cell_id_list.append(adata_group_prot.obs.cell_id)
# # print(sample_id_list)
# # print(len(adata_group_prot.obs))
# # print(len(sample_id_list))

# for file in group:
#     print("Current sample: ", file)
#     #adata = anndata.read("/panfs/biopan04/unr_thailand/kb_out/"+prefix+'/counts_unfiltered/adata.h5ad')

#     if file == first_sample:
#         print("First Sample")

#     else:
#         adata_p = anndata.read("/panfs/biopan04/unr_thailand/cell_ranger_CSP_h5ad_new_antibody_list/"+file)

#         sample_id_list.append([file.split(".")[0]] * len(adata_p.obs))
#         cell_id_list.append(adata_p.obs.cell_id)

#         adata_group_prot = ad.concat([adata_group_prot, adata_p])
#         print(adata_group_prot)

# adata_group_prot.var["gene_id"] = gene_id
# print("")
# print("adata_group_prot")
# print(adata_group_prot)
# print(adata_group_prot.obs["cell_id"])

# #Adds sample ids to protein adata
# sample_id_list = flatten(sample_id_list)
# adata_group_prot.obs["sample_id"] = sample_id_list
# cell_id_list = flatten(cell_id_list)
# print(len(cell_id_list))
# print(cell_id_list[0])
# adata_group_prot.obs["cell_id_2"] = cell_id_list

# print(adata_group_prot.obs.cell_id_2)
# print(len(adata_group_prot.obs.cell_id_2))

# print(adata_group_prot.obs.sample_id)
# print(len(adata_group_prot.obs.sample_id))

# print(adata_group_prot.obs_names)
# print(len(adata_group_prot.obs_names))

# print(adata_group_prot[["0", "1"],[]])




# #fixes cell id names
# cell_id = adata_group_prot.obs
# print(cell_id)
# cell_id = cell_id["cell_id"].to_list()
# cell_id = [k[:-2] for k in cell_id] #remove -1 off the cell id so it matches the gene format
# adata_group_prot.obs["cell_id"] = cell_id

# adata_prot = adata_group_prot
# print("adata_prot")
# print(adata_prot)
# print(adata_prot.obs["cell_id"])
# print(adata_prot.obs.cell_id)

print("")
print("gene data")
print(adata)
print(adata.obs_names)
# print(adata.obs.condition)



































### Subset by treatment group if needed ####
if outcome_group not in "All":
    print("Selecting subgroup....")
    adata = adata[adata.obs.condition == outcome_group]
    print(adata)



### Recluster Leiden Clustering (if needed) ###
#Default clustering is 0.5, can recluster with a different number here
# scanpy.pp.neighbors(adata, use_rep = 'X_scVI')
# scanpy.pp.neighbors(adata, n_neighbors=30, n_pcs=20)
# scanpy.tl.umap(adata)
scanpy.tl.leiden(adata, resolution = 0.4)



### Plot Leiden Clustering ###
with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    scanpy.pl.umap(adata, color='leiden', ax=ax)
    plt.savefig(out_dir+"UMAP_leiden_"+outcome_group+".png")



### Calculates Marker Genes ###
#subset by leiden group
groups = adata.obs["leiden"]
groups = set(groups.to_list())
groups = list(groups)
groups.sort(key=int)
print(groups)


#Use t-test to find differentailly expressed genes for each group
#Gene 
scanpy.tl.rank_genes_groups(adata, 'leiden', method="t-test")
markers = scanpy.get.rank_genes_groups_df(adata, groups)
markers_gene = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
print(markers_gene)

with plt.rc_context():
    scanpy.pl.rank_genes_groups(adata)
    plt.savefig(out_dir+"rank_genes_groups_gene_res_1.0_"+outcome_group+".png")

# #Protein
# #Artifically map gene leiden groups to protein data
# artifical_leiden = []
# adata_list = []
# valid_prot_group = [] #Holds the list of groups 
# for group_num in groups:
#     print("Analyzing group: ", group_num)
#     adata_subset = adata[adata.obs["leiden"] == group_num]
#     group_total_prot = 0

#     #Subsets futher by sample id to aviod id reuse across different samples
#     # print(type(adata_subset.obs.sample_id))

#     for sample_id in set(adata_subset.obs.sample_id):
#         print("On sample id... ", sample_id)
#         adata_subset_sample = adata_subset[adata_subset.obs["sample_id"] == sample_id]
#         print(adata_prot.obs.cell_id)
#         print(adata_prot)
#         print(adata_prot.obs.sample_id)
#         # prot_id_list = adata_prot.obs.sample_id
#         # print(prot_id_list)
#         # with open(out_dir+group_num+"Protien_df_sample_id_debug.csv", 'w') as f:
#         #     for line in prot_id_list:
#         #         f.write("%s\n" % line)

#         adata_prot_sample = adata_prot[adata_prot.obs["sample_id"] == sample_id]
#         print(adata_prot_sample)
#         print(adata_prot_sample.obs.cell_id)
#         # print("Subsets")
#         # print(adata_subset_sample)
#         # print(adata_prot_sample.obs)
#         # prot_df = adata_prot_sample.to_df()
#         # print("PROT LIST")
#         # prot_id_list = adata_prot_sample.obs["cell_id_2"]
#         # print(prot_id_list)
#         # with open(out_dir+group_num+"Protien_df_debug.csv", 'w') as f:
#         #     for line in prot_id_list:
#         #         f.write("%s\n" % line)
#         # prot_df.set_index("cell_id")
#         # print(prot_df.head())
#         # prot_df.to_csv(out_dir+group_num+"Protien_df_debug.csv", sep='\t')

#         # gene_df = adata_subset.to_df()
#         # # gene_df["cell_id"] = adata_subset.obs["cell_id"]
#         # # gene_df.set_index("cell_id")
#         # print(gene_df.head())
#         # gene_df.to_csv(out_dir+group_num+"Gene_df_debug.csv", sep='\t')

#         # print("GENE LIST")
#         # gene_id_list = adata_subset_sample.obs_names
#         # print(gene_id_list)
#         # with open(out_dir+group_num+"Gene_df_debug.csv", 'w') as f:
#         #     for line in gene_id_list:
#         #         f.write("%s\n" % line)

#         print(adata_prot_sample.obs.cell_id)

#         #find cell ids that have protein data and are in the subgroup
#         matching_ids = []
#         for item in adata_subset_sample.obs_names:
#             if item in adata_prot_sample.obs.cell_id.values:
#                 matching_ids.append(item)

#         # print(matching_ids)
#         # print("matching ids")
#         # print(len(matching_ids))
#         group_total_prot = len(matching_ids) + group_total_prot
#         #Everything is in the list twice idk why but this fixes it
#         # matching_ids = list(set(matching_ids))


#         #Protein data subset for leiden group
#         id_match_adata_prot = adata_prot_sample[adata_prot_sample.obs.cell_id.isin(matching_ids)]
#         # print("id_match_adata_prot")
#         # print(id_match_adata_prot)
        
#         #adds a list of the group number with length obs to the master list artifical_leiden
#         artifical_leiden.append([group_num] * len(id_match_adata_prot.obs))

#         #adds anndata to list
#         adata_list.append(id_match_adata_prot)
#     print("Total matching protein for group is ", group_total_prot)

#     #Prevents issues with rank gene groups when there is not > 1 cell ids that match
#     if group_total_prot > 1:
#         print("Valid Number of Protien ids found ", group_total_prot)
#         valid_prot_group.append(group_num)

# print("Valid protein groups list ")
# print(valid_prot_group)

# print("adata prot below")
# print(adata_prot)
# print(adata_prot.obs.cell_id)

#combines all anndata group obj
# adata_prot = ad.concat(adata_list)

# artifical_leiden = flatten(artifical_leiden)
# adata_prot.obs["leiden"] = artifical_leiden
# adata_prot.var["gene_id"] = gene_id

# print("Protein adata")
# print(adata_prot)
# print(adata_prot.obs)
# print(adata_prot.var)

#Data needs to be normalized for rank_gene_groups
# scanpy.pp.normalize_total(adata_prot, target_sum = 1e4)
# scanpy.pp.log1p(adata_prot)

# #Identifies a list of marker genes for each protein group
# scanpy.tl.rank_genes_groups(adata_prot, 'leiden', groups=valid_prot_group, method="t-test")
# markers = scanpy.get.rank_genes_groups_df(adata_prot, valid_prot_group)
# markers_prot = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
# print(markers_prot)

#Adds a column to the markers df with the cell type
# markers_prot["Cell_Type"] = ""
# for index, row in markers_prot.iterrows():
#     markers_prot.loc[index, "Cell_Type"] = index_gene_id_key.get(markers_prot.loc[index, "names"])

# print(markers_prot)

# #Identifies a list of marker genes for the lowly expressed markers
# markers_prot_low = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges < -.5)]
# print(markers_prot_low)

# #Adds cell type to markers df
# markers_prot_low["Cell_Type"] = ""
# for index, row in markers_prot_low.iterrows():
#     markers_prot_low.loc[index, "Cell_Type"] = index_gene_id_key.get(markers_prot_low.loc[index, "names"])

# print(markers_prot_low)


# with plt.rc_context():
#     scanpy.pl.rank_genes_groups(adata_prot)
#     plt.savefig(out_dir+"rank_genes_groups_prot_res_1.0_"+outcome_group+".png")





### Calculate Automated Cell Type ####

calculated_cell_types_gene = []
calculated_cell_types_prot_high = []
calculated_cell_types_prot_low = []
calculated_cell_types_gene_score = []
calculated_cell_types_prot_high_score = []
calculated_cell_types_prot_low_score = []
calculated_cell_types_consensus = []

for group_num in groups:
    print("Group number: ", group_num)
    adata_subset = adata[adata.obs["leiden"] == group_num]
    # id_match_adata_prot = adata_prot[adata_prot.obs["leiden"] == group_num]

    #Creates a list of the gene ids
    # var_list = id_match_adata_prot.var.gene_id
    # var_list = var_list.to_dict()
    # var_list = var_list.values()

    # #Sums the cell surface protein data and creates a df
    # id_match_adata_prot = scanpy.AnnData(X = id_match_adata_prot.X.sum(axis = 0), var = var_list)
    # df = id_match_adata_prot.to_df()
    # df = df.transpose()
    # df["Id"] = list(var_list) #Adds gene ids
    # df = df.loc[(df!=0).any(axis=1)]
    # df = df.sort_values(by=['0'], ascending=False)

 
    # #Find the protein based cell types that match the group
    # df["Cell_Type"] = ""
    # for index, row in df.iterrows():
    #     if row[1] in cell_types_prot_pos:
    #         df.at[index, 'Cell_Type'] = cell_types_prot_pos.get(row[1])
    #     else:
    #         df.at[index, 'Cell_Type'] = ['NA']
   
   
#    #Checks if the df is empty and sets an NA value
#     if len(df) == 0:
#         df = pd.DataFrame([["No Proteins", 0, np.array(['NA'])]], columns=["Prot", "0", "Cell_Type"])
#         df.index = df["Prot"]
#         df = df.drop('Prot', axis=1)


    #creates df with sum of the counts for each gene
    count_list = adata_subset.X.sum(axis = 0).tolist()
    count_list = [item for items in count_list for item in items]

    dict = {"Gene": adata_subset.var_names.to_list(), "Sum_Counts": count_list} 
    df_gene = pd.DataFrame(dict)
    df_gene.index = df_gene["Gene"]
    df_gene = df_gene.sort_values(by=['Sum_Counts'], ascending=False)

    #Gene set cell type match
    for item in adata_subset.var_names:
        if item in cell_types_gene:
            df_gene.at[item, 'Cell_Type'] = cell_types_gene[item]
        else:
            df_gene.at[item, 'Cell_Type'] = ['NA']
    
    # Created file with a list of all the genes found in the data (includes NA)
    # df.to_csv(out_dir+group_num+"_top_CSP_cell_type_res_1.0_"+outcome_group+".csv")
    df_gene.to_csv(out_dir+group_num+"_top_GEX_cell_type_res_1.0_"+outcome_group+".csv")


    #Counts number of occurances of the marker genes for the cell types in the data
    #Genes
    cell_type_gene_occurrences = {} #hold the total number of times a marker gene occurances in the data for each cell type
    for index, row in df_gene.iterrows():
        if row[1] > 0: #Only counts genes that are present in the dataset with counts > 0
            for item in row[2]:
                if item in cell_type_gene_occurrences:
                    cell_type_gene_occurrences[item] = cell_type_gene_occurrences.get(item)+1
                else:
                    cell_type_gene_occurrences[item] = 1
    del cell_type_gene_occurrences["NA"]

    # #Pos proteins
    # cell_type_protein_occurrences = {} #hold the total number of times a marker gene occurances in the data for each cell type
    # for index, row in df.iterrows():
    #     if row[0] > 0:
    #         for item in row[2]:
    #             if item in cell_type_protein_occurrences:
    #                 cell_type_protein_occurrences[item] = cell_type_protein_occurrences.get(item)+1
    #             else:
    #                 cell_type_protein_occurrences[item] = 1
    # if bool(cell_type_protein_occurrences) == True: #checks it is not empty before deleting to prevent errors
    #     del cell_type_protein_occurrences["NA"]
    
    
    # Created file with a list of all the genes with cell types found in the data
    gene_df_filtered = df_gene[~df_gene.Cell_Type.str.contains("NA", regex=False)]
    # prot_df_filtered = df[~df.Cell_Type.str.contains("NA", regex=False)]
    gene_df_filtered.to_csv(out_dir+group_num+"_top_GEX_cell_type_no_NA_res_1.0.csv")
    # prot_df_filtered.to_csv(out_dir+group_num+"_top_CSP_cell_type_no_NA_res_1.0.csv")
  


    #Calcuate num of genes ("points") per cell type
    #Protein highly expressed
    # cell_type_count_prot_high = {}
    # markers_prot_group = markers_prot[(markers_prot['group'] == group_num)]  #gets list of markers that belong to that leiden group
   
    # for index, row in markers_prot_group.iterrows():
    #     if row[6] in cell_types_prot_pos:
    #         for item in cell_types_prot_pos[row[6]]: #loops through the list of cell types for that marker
    #             if item in cell_type_count_prot_high: #if the cell type already exists in the dictionary it +1
    #                 cell_type_count_prot_high[item] = cell_type_count_prot_high.get(item)+1
    #             else: #add cell type to dictionary if it does not exist
    #                 cell_type_count_prot_high[item] = 1

    # #Protein lowly expressed
    # cell_type_count_prot_low = {}
    # markers_prot_low_group = markers_prot_low[(markers_prot_low['group'] == group_num)]    
    # for index, row in markers_prot_low_group.iterrows():
    #     if row[6] in cell_types_prot_neg:
    #         for item in cell_types_prot_neg[row[6]]:
    #             if item in cell_type_count_prot_low:
    #                 cell_type_count_prot_low[item] = cell_type_count_prot_low.get(item)+1
    #             else:
    #                 cell_type_count_prot_low[item] = 1
    
    #Gene   
    cell_type_count_gene = {}
    markers = markers_gene[(markers_gene['group'] == group_num)]

    for index, row in markers.iterrows():
        if row[1] in cell_types_gene:
            for item in cell_types_gene[row[1]]:
                if item in cell_type_count_gene:
                    cell_type_count_gene[item] = cell_type_count_gene.get(item)+1
                else:
                    cell_type_count_gene[item] = 1




    #Calculate Score
    #number of time the cell type occurs in the marker genes divided by the number of times it occurs in the dataset group
    print("Gene marker gene count")
    print(cell_type_count_gene)
    print("Gene occurrences")
    print(cell_type_gene_occurrences)

    cell_type_gene_score = {}
    for cell_type, occurrences in cell_type_gene_occurrences.items():
        if cell_type_count_gene.get(cell_type): #key exists
            cell_type_gene_score[cell_type] = int(cell_type_count_gene[cell_type]) / int(occurrences)
            # print("Score is ", cell_type_gene_score[cell_type], " for ", str(cell_type_count_gene[cell_type]), " over ", occurrences)
    print("Gene score")
    print(cell_type_gene_score)

    # print("Protein high marker gene count")
    # print(cell_type_count_prot_high)
    # print("Protein high occurences")
    # print(cell_type_protein_occurrences)
    # cell_type_prot_high_score = {}
    # for cell_type, occurrences in cell_type_protein_occurrences.items():
    #     if cell_type_count_prot_high.get(cell_type): #key exists
    #         cell_type_prot_high_score[cell_type] = int(cell_type_count_prot_high[cell_type]) / int(occurrences)
    #         # print("Score is ", cell_type_prot_high_score[cell_type], " for ", str(cell_type_count_prot_high[cell_type]), " over ", occurrences)
    # print("Protein high score")   
    # print(cell_type_prot_high_score)


    # print("Protein low marker gene count")
    # print(cell_type_count_prot_low)
    # print("Protein low occurences")
    # print(cell_type_protein_occurrences)

    # cell_type_prot_low_score = {}
    # for cell_type, occurrences in cell_type_protein_occurrences.items():
    #     if cell_type_count_prot_low.get(cell_type): #key exists
    #         cell_type_prot_low_score[cell_type] = int(cell_type_count_prot_low[cell_type]) / int(occurrences)
    #         # print("Score is ", cell_type_prot_low_score[cell_type], " for ", str(cell_type_count_prot_low[cell_type]), " over ", occurrences)
    # print("Protein low score")   
    # print(cell_type_prot_low_score)

            


    #Find cell type with the highest score
    #Protein
    # if bool(cell_type_prot_high_score) == False: #if dict cell_type_prot_high_score is empty
    #     max_val_prot_high = "NA_"+group_num
    #     max_val_prot_high_score = 0
    # else:
    #     score_vals = find_maxes(cell_type_prot_high_score)
    #     max_val_prot_high = score_vals[0]
    #     max_val_prot_high_score = score_vals[1]
    #     if max_val_prot_high == "NA" and len(cell_type_prot_high_score) > 1:
    #         #find second best because even though NA has the most points another real cell type is present
    #         del cell_type_prot_high_score['NA']
    #         score_vals = find_maxes(cell_type_prot_high_score)
    #         max_val_prot_high = score_vals[0]
    #         max_val_gene_prot_high = score_vals[1]
    #     elif max_val_prot_high == "NA":
    #         max_val_prot_high = max_val_prot_high+"_"+group_num
    #     else: #only NAs for cell type
    #         max_val_prot_high = max_val_prot_high
    # print("max val prot high ", max_val_prot_high)


    # if bool(cell_type_prot_low_score) == False:
    #     max_val_prot_low = "NA_"+group_num
    #     max_val_prot_low_score = 0
    # else:
    #     score_vals = find_maxes(cell_type_prot_low_score)
    #     max_val_prot_low = score_vals[0]
    #     max_val_prot_low_score = score_vals[1]       
    #     if max_val_prot_low == "NA" and len(cell_type_prot_low_score) > 1:
    #         #find second best because even though NA has the most points another real cell type is present
    #         del cell_type_prot_low_score['NA']
    #         score_vals = find_maxes(cell_type_prot_low_score)
    #         max_val_prot_low = score_vals[0]
    #         max_val_gene_prot_low = score_vals[1]
    #     elif max_val_prot_low == "NA":
    #         max_val_prot_low = max_val_prot_low+"_"+group_num
    #     else: #only NAs for cell type
    #         max_val_prot_low = max_val_prot_low
    # print("max val prot low ", max_val_prot_low)

    #Gene
    if bool(cell_type_gene_score) == False:
        max_val_gene = "NA_"+group_num
        max_val_gene_score = 0
    else:
        score_vals = find_maxes(cell_type_gene_score)
        max_val_gene = score_vals[0]
        max_val_gene_score = score_vals[1]
        if max_val_gene == "NA" and len(cell_type_gene_score) > 1:
            #find second best because even though NA has the most points another real cell type is present
            del cell_type_gene_score['NA']
            score_vals = find_maxes(cell_type_gene_score)
            max_val_gene = score_vals[0]
            max_val_gene_score = score_vals[1]
        elif max_val_gene == "NA": #only NAs for cell type
            max_val_gene = max_val_gene+"_"+group_num
        else:
            max_val_gene = max_val_gene
    print("max val gene ", max_val_gene)

    #makes sure group names are unique
    # #Protein
    # if max_val_prot_high in calculated_cell_types_prot_high:
    #     max_val_prot_high = max_val_prot_high+" "+group_num

    # if max_val_prot_low in calculated_cell_types_prot_low:
    #     max_val_prot_low = max_val_prot_low+" "+group_num
    
    # #Gene
    # if max_val_gene in calculated_cell_types_gene:
    #     max_val_gene = max_val_gene+" "+group_num

    # calculated_cell_types_prot_high.append(max_val_prot_high)
    # calculated_cell_types_prot_high_score.append(max_val_prot_high_score)
    # calculated_cell_types_prot_low.append(max_val_prot_low)
    # calculated_cell_types_prot_low_score.append(max_val_prot_low_score)
    calculated_cell_types_gene.append(max_val_gene)
    calculated_cell_types_gene_score.append(max_val_gene_score)


# print("Protein high calculated types")
# print(calculated_cell_types_prot_high)
# print(calculated_cell_types_prot_high_score)
# print("Protein low calculated types")
# print(calculated_cell_types_prot_low)
# print(calculated_cell_types_prot_low_score)
print("Gene calculated types")
print(calculated_cell_types_gene)
print(calculated_cell_types_gene_score)

#Allows for you to only consider the gene expression data
if only_gene_expr is True:

    for i,j in zip(calculated_cell_types_gene,calculated_cell_types_gene_score):
        calculated_cell_types_consensus.append(str(i)+" ("+str(j)+")")

    print("Gene Expression Only Results Below:")
    print(calculated_cell_types_consensus)


else:
    ### Resolves discrepancies between the different cell types ####
    print("Comparison below: ")
    for i in range(0, len(calculated_cell_types_prot_high)):
        print(" ")
        print(calculated_cell_types_gene[i], " ", calculated_cell_types_prot_high[i], " ", calculated_cell_types_prot_low[i])
        cell_type_options = [calculated_cell_types_gene[i],calculated_cell_types_prot_high[i], calculated_cell_types_prot_low[i]]
        #Agree
        if calculated_cell_types_gene[i] == calculated_cell_types_prot_high[i] == calculated_cell_types_prot_low[i]:
            score = find_score(calculated_cell_types_prot_high[i], i)
            calculated_cell_types_consensus.append(calculated_cell_types_prot_high[i]+" ("+str(score)+")")
            print("All cell types were the same, appending ", calculated_cell_types_prot_high[i])

        #1 of the options is NA
        elif sum('NA_' in s for s in cell_type_options) == 1:
            #Other cell types match
            if len(set(cell_type_options)) == 2:
                list_types = [x for x in cell_type_options if 'NA_' not in x]
                # calculated_cell_types_consensus.append(list_types[0])
                score = find_score(list_types[0], i)
                calculated_cell_types_consensus.append(list_types[0]+" ("+str(score)+")")
                print("1 NA in list and others match appending ",list_types[0])
            #Other cell types do not match, returns type with a larger score
            else:
                #check if any of the cell types are 2 appended together (ex "MC | Eos")
                if any("|" in s for s in cell_type_options):
                    print("One type contains multible options")
                    index = 0
                    cell_type_dict_score = {}
                    for item in cell_type_options:
                        if "NA_" not in item:
                            if index == 0:
                                cell_type_dict_score[item] = calculated_cell_types_gene_score[i]
                            elif index == 1:
                                cell_type_dict_score[item] = calculated_cell_types_prot_high_score[i]
                            elif index == 2:
                                cell_type_dict_score[item] = calculated_cell_types_prot_low_score[i]
                        index+=1

                    cell_type_dict = {}
                    for key, value in cell_type_dict_score.items():
                        if "|" in key: #handles cases with multible top cell types
                            for subitem in  key.split("|"):
                                if subitem in cell_type_dict:
                                    cell_type_dict[subitem] = (float(cell_type_dict.get(subitem))+float(value))
                                else:
                                    cell_type_dict[subitem] = value
                        else:
                            if key in cell_type_dict:
                                cell_type_dict[key] = (float(cell_type_dict.get(key))+float(value))
                            else:
                                cell_type_dict[key] = value
                    
                    print(cell_type_dict)
                    score_vals = find_maxes(cell_type_dict)
                    max_val_gene = score_vals[0]
                    score = score_vals[1]
                    print("After dividing multi and adding scores, chosing ", max_val_gene, " with a score of ", score)
                    calculated_cell_types_consensus.append(max_val_gene+" ("+str(score)+")") 

                else:
                    index = 0
                    cell_type_dict = {}
                    for item in cell_type_options:
                        if "NA_" not in item:
                            if index == 0:
                                cell_type_dict[item] = calculated_cell_types_gene_score[i]
                            elif index == 1:
                                cell_type_dict[item] = calculated_cell_types_prot_high_score[i]
                            elif index == 2:
                                cell_type_dict[item] = calculated_cell_types_prot_low_score[i]
                        index+=1

                    score_vals = find_maxes(cell_type_dict)
                    max_val_gene = score_vals[0]
                    score = score_vals[1]
                    calculated_cell_types_consensus.append(max_val_gene+" ("+str(round(score, 2))+")")
                    print("1 NA, 2 other types do not match, appending ", max_val_gene, " with a score of ", score)           

        #2 of the options is NA
        #elif cell_type_options.count('NA_') == 2:
        elif sum('NA_' in s for s in cell_type_options) == 2:
            #handles cases where the non NA option has multible top scores
            #check if any of the cell types are 2 appended together (ex "MC | Eos")
            if any("|" in s for s in cell_type_options):
                print("One type contains multible options")
                index = 0
                cell_type_dict_score = {}
                for item in cell_type_options:
                    if "NA_" not in item:
                        if index == 0:
                            cell_type_dict_score[item] = calculated_cell_types_gene_score[i]
                        elif index == 1:
                            cell_type_dict_score[item] = calculated_cell_types_prot_high_score[i]
                        elif index == 2:
                            cell_type_dict_score[item] = calculated_cell_types_prot_low_score[i]
                    index+=1
                cell_type_dict = {}
                for key, value in cell_type_dict_score.items():
                    if "|" in key: #handles cases with multible top cell types
                        # print("its a multi ", key)
                        for subitem in  key.split("|"):
                            # print("subitem: ", subitem)
                            if subitem in cell_type_dict:
                                cell_type_dict[subitem] = (float(cell_type_dict.get(subitem))+float(value))
                            else:
                                cell_type_dict[subitem] = value
                    else:
                        if key in cell_type_dict:
                            cell_type_dict[key] = (float(cell_type_dict.get(key))+float(value))
                        else:
                            cell_type_dict[key] = value
                
                print(cell_type_dict)
                score_vals = find_maxes(cell_type_dict)
                max_val_gene = score_vals[0]
                score = score_vals[1]
                print("After dividing multi and adding scores, chosing ", max_val_gene, " with a score of ", score)
                calculated_cell_types_consensus.append(max_val_gene+" ("+str(round(score, 2))+")")
            else:
                list_types = [x for x in cell_type_options if 'NA_' not in x]
                score = find_score(list_types[0], i)
                calculated_cell_types_consensus.append(list_types[0]+" ("+str(round(score, 2))+")")
                print("2 options is NA, 1 is NOT  ", list_types[0])

        #3 cell type options that do not all agree, 2 types match
        elif len(set(cell_type_options)) == 2:
            occurence_count = Counter(cell_type_options)
            most_common_cell_type  = occurence_count.most_common(1)[0][0]
            score = find_score(most_common_cell_type, i)

            print("3 valid cell types, 2 agree chosing ", most_common_cell_type, " with a score of ", score)
            calculated_cell_types_consensus.append(most_common_cell_type+ " ("+str(round(score, 2))+")")

        #3 cell type options that do not all agree, none types match
        elif len(set(cell_type_options)) == 3:
            #handles cases where the non NA option has multible top scores
            #check if any of the cell types are 2 appended together (ex "MC | Eos")
            if any("|" in s for s in cell_type_options):
                print("One type contains multible options")
                index = 0
                cell_type_dict_score = {}
                for item in cell_type_options:
                    if "NA_" not in item:
                        if index == 0:
                            cell_type_dict_score[item] = calculated_cell_types_gene_score[i]
                        elif index == 1:
                            cell_type_dict_score[item] = calculated_cell_types_prot_high_score[i]
                        elif index == 2:
                            cell_type_dict_score[item] = calculated_cell_types_prot_low_score[i]
                    index+=1

                cell_type_dict = {}
                for key, value in cell_type_dict_score.items():
                    if "|" in key: #handles cases with multible top cell types
                        for subitem in  key.split("|"):
                            if subitem in cell_type_dict:
                                cell_type_dict[subitem] = (float(cell_type_dict.get(subitem))+float(value))
                            else:
                                cell_type_dict[subitem] = value
                    else:
                        if key in cell_type_dict:
                            cell_type_dict[key] = (float(cell_type_dict.get(key))+float(value))
                        else:
                            cell_type_dict[key] = value
                
                print(cell_type_dict)
                score_vals = find_maxes(cell_type_dict)
                max_val_gene = score_vals[0]
                score = score_vals[1]
                print("3 valid options with multi and adding scores, chosing ", max_val_gene, " with a score of ", score)
                calculated_cell_types_consensus.append(max_val_gene+" ("+str(round(score, 2))+")") 
            else:
                index = 0
                cell_type_dict = {}
                for item in cell_type_options:
                    if index == 0:
                        cell_type_dict[item] = calculated_cell_types_gene_score[i]
                    elif index == 1:
                        cell_type_dict[item] = calculated_cell_types_prot_high_score[i]
                    elif index == 2:
                        cell_type_dict[item] = calculated_cell_types_prot_low_score[i]
                    index+=1

                score_vals = find_maxes(cell_type_dict)
                max_val_gene = score_vals[0]
                score = score_vals[1]
                print("All 3 different valid cell types, chosing ", max_val_gene, " with a score of ", score)
                calculated_cell_types_consensus.append(max_val_gene+" ("+str(round(score, 2))+")")
            
        #Catch all, in case something slips through
        else:
            calculated_cell_types_consensus.append(str(calculated_cell_types_gene[i])+" or "+str(calculated_cell_types_prot_high[i])+" or "+str(calculated_cell_types_prot_low[i]) )
            print("in the catch all")
            
# #resolves prot and gene
# for i in range(0, len(calculated_cell_types_prot_high)):
#     print(calculated_cell_types_gene[i], " ", calculated_cell_types_prot_high[i], " ", calculated_cell_types_prot_low)
#     #Agree
#     if calculated_cell_types_gene[i] == calculated_cell_types_prot[i]:
#         calculated_cell_types_consensus.append(calculated_cell_types_prot[i])
#     #Gene is NA, Prot is resolved
#     elif ('NA_' in calculated_cell_types_gene[i]) and ('NA_' not in calculated_cell_types_prot[i]):
#         calculated_cell_types_consensus.append(calculated_cell_types_prot[i])
#     #Gene is resolved, Prot is NA
#     elif ('NA_' not in calculated_cell_types_gene[i]) and ('NA_' in calculated_cell_types_prot[i]):
#         calculated_cell_types_consensus.append(calculated_cell_types_gene[i])
#     #both not NAs but disagree   
#     else:
#         calculated_cell_types_consensus.append(str(calculated_cell_types_gene[i])+" or "+str(calculated_cell_types_prot[i]))



#### Rename groups with calculated names #### 

# only_name_calculated_cell_types_consensus = []
unique_calculated_cell_types_consensus = []
count = 0
for cell_type in calculated_cell_types_consensus:
    # print(cell_type.split(" (")[0])
    if cell_type in unique_calculated_cell_types_consensus:
        # print(cell_type+" is already in ")
        unique_calculated_cell_types_consensus.append(cell_type+" "+str(count))
    else:
        unique_calculated_cell_types_consensus.append(cell_type)
        # unique_calculated_cell_types_consensus.append(cell_type.split(" (")[0])


    count+=1

print("Orignial Names")
print(calculated_cell_types_consensus)
# print(only_name_calculated_cell_types_consensus)
print("Unique Names Below")
print(unique_calculated_cell_types_consensus)

calculated_cell_types_consensus = unique_calculated_cell_types_consensus


#Renames groups
print("Final cell types with scores")
print(calculated_cell_types_consensus)
adata.rename_categories('leiden', calculated_cell_types_consensus)

#Saves new h5ad file with updated names 
adata.write_h5ad(out_dir+"adata_w_leiden_groups_res_1.0_"+outcome_group+"_GENE_ONLY.h5ad")

#Replots the leiden clustering with new groups
with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    scanpy.pl.umap(adata, color='leiden', ax=ax)
    plt.savefig(out_dir+"UMAP_leiden_automated_names_res_1.0_"+outcome_group+"_GENE_ONLY.png", bbox_inches='tight')
with plt.rc_context():
    fig, ax = plt.subplots(figsize=(10, 7))
    scanpy.pl.umap(adata, color='leiden', legend_loc='on data', frameon=False, legend_fontsize=10, legend_fontoutline=2)
    plt.savefig(out_dir+"UMAP_leiden_automated_names_labeled_on_plot_res_1.0_"+outcome_group+"_GENE_ONLY.png", bbox_inches='tight')





#### Dot Plots ####

# #creates a dot plot with all the marker genes from the marker gene file from the colaborators
# gene_list = list(cell_types_gene.keys())

# marker_gene_file = [] #hold list of marker genes currently in this file
# for marker_gene in gene_list:
#     if marker_gene in adata.var_names.to_list():
#         marker_gene_file.append(marker_gene)
# # print(marker_gene_file)

# with plt.rc_context():
#     # fig, ax = plt.subplots(figsize=(10, 7))
#     scanpy.pl.dotplot(adata, marker_gene_file, groupby="leiden")
#     plt.savefig(out_dir+"dot_plot_marker_genes.png")


# #dot plot with genes that the gene ttest found and that were in the cell type marker gene list from the file
# marker_gene_list_ttest = markers['names'].tolist()
# intersecting_genes = list((Counter(marker_gene_file) & Counter(marker_gene_list_ttest)).elements())

# with plt.rc_context():
#     # fig, ax = plt.subplots(figsize=(10, 7))
#     scanpy.pl.dotplot(adata, intersecting_genes, groupby="leiden")
#     plt.savefig(out_dir+"dot_plot_marker_genes_from_ttest_also_marker.png")

