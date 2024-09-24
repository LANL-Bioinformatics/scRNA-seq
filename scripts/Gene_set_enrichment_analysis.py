# Runs gene set enrichment on results from the differential gene expression

import numpy as np
import matplotlib.pyplot as plt
import argparse
import warnings
warnings.filterwarnings("ignore")  # plotnine has a lot of MatplotlibDeprecationWarning's
import pandas as pd
import gseapy as gp
import os


# Define the parser
parser = argparse.ArgumentParser(description='Cell Composition Analysis')

parser.add_argument('--n_top_term_enrich_plot', dest='n_top_term_enrich_plot', type=int, default=5)
parser.add_argument('--dotplot_cutoff', dest='dotplot_cutoff', type=float, default=0.25)
parser.add_argument('--enrich_terms', dest='enrich_terms', nargs='*', default= ["GO_Biological_Process_2023","GO_Molecular_Function_2023"])
parser.add_argument('--output_folder', action="store", dest='out_dir', required=True)

# Parse the command line arguments and store in args
args = parser.parse_args()

print("*** Code parameters ***")
print("Number of top terms displayed in enrichment plots: ", args.n_top_term_enrich_plot)
print("Dotplot cutoff: ", args.dotplot_cutoff)
print("Gene sets for enrichment: ", args.enrich_terms)
print("Output folder: ", args.out_dir)
print("***********************")

out_dir = args.out_dir

# Gets list of files with the right extenstion
dir_list = os.listdir(out_dir)
files = [f for f in dir_list if "_all_genes_differential_gene_expr_results.csv" in f]
files.sort()

for file in files:
    
    print("\nProcessing file: ", file)
    prefix = file.split("_all_genes_differential_gene_expr_results.csv")[0]


    # Load data
    print("\nLoading data and calculating test statistic... ")
    df = pd.read_csv(out_dir+file)

    # Data setup and calculate statistic
    df.rename(columns={'Unnamed: 0':'Gene'}, inplace=True ) 

    df['Rank'] = -np.log10(df.padj)*df.log2FoldChange

    df = df.sort_values('Rank', ascending = False).reset_index(drop = True)

    ranking = df[['Gene', 'Rank']].dropna().sort_values('Rank', ascending=False)

    # print(gp.get_library_name())
    for term in args.enrich_terms:
        print("\nAnalyzing with term: ", term)
        print("\nRunning GSEA...")
        pre_res = gp.prerank(rnk = ranking,
                            gene_sets = term,
                            seed = 6, 
                            permutation_num = 100) #INCRESE TO AT LEAST 100


        try:
            ax = gp.plot.dotplot(pre_res.res2d,
                        column="FDR q-val",
                        title=term,
                        cmap=plt.cm.viridis,
                        size=6, # adjust dot size
                        figsize=(4,5), 
                        cutoff=args.dotplot_cutoff, 
                        show_ring=False, 
                        ofname=out_dir+prefix+"_"+term+"_dotplot_enrichment.png")
        except ValueError as e: #moves on if nothing meets the cutoff
            print(e)
            print("Skipping dotplot")

        terms = pre_res.res2d.Term
        axs = pre_res.plot(terms=terms[1], ofname=out_dir+prefix+"_"+term+"_top_term_enrichment.png")

        axs = pre_res.plot(terms=terms[1:args.n_top_term_enrich_plot+1],
                        #legend_kws={'loc': (1.2, 0)}, # set the legend loc
                        show_ranking=True, # whether to show the second yaxis
                        figsize=(3,4),
                        ofname=out_dir+prefix+"_"+term+"_top_"+str(args.n_top_term_enrich_plot)+"_term_enrichment.png"
                        )
        print("\nTop ", str(args.n_top_term_enrich_plot), " results below... Full results located in ", prefix+"_"+term+"_enrichment_scores.csv")
        print(pre_res.res2d.head(args.n_top_term_enrich_plot))
        pre_res.res2d.to_csv(out_dir+prefix+"_"+term+"_enrichment_scores.csv")



