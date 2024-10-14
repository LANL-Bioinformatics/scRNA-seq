## scRNA-seq

## Introduction

## Workflow

This workflow is automated using a Nextflow workflow described in `sc_pipeline.nf`.

### Requirements

Running `sc_pipeline.nf` requires an environment with Nextflow installed, in addition to container management software such as Singularity or Charliecloud (see [Containers](#containers)). This workflow was developed and tested using Nextflow v24.04.4.


### Containers

All necessary scripts, third-party software, and dependencies are included in a Docker image hosted at https://hub.docker.com/repository/docker/apwat/10x_sc/. The Dockerfile and Conda environment file used to create the image are archived in this repository. Nextflow will attempt to pull and convert this image to one of the runtimes below, and then run its processes in containers.

#### Singularity

The default container software used to run processes in `sc_pipeline.nf`. 
Setting the environment variable `NXF_SINGULARITY_CACHEDIR` will control where images are downloaded.

```
Tested with SingularityCE 3.9.5 and up -- older versions of Singularity may be unable to properly set up containers.
```

#### Charliecloud

*Use not recommended until https://github.com/nextflow-io/nextflow/pull/5300 is included in a Nextflow release.*

To instead use Charliecloud as a container runtime, comment out the `singularity` scope in `nextflow.config`, and uncomment the `charliecloud` scope.

As of Nextflow v24.04.4, workflows with multiple processes running in the same Charliecloud container will attempt multiple pulls of the same image, which fail. To resolve this, download the image ahead of time:

```
export CH_IMAGE_STORAGE=/path/to/image_storage/
export NXF_CHARLIECLOUD_CACHEDIR=/path/to/image_storage/
ch-image pull apwat/10x_sc:0.6
```
-and then run the workflow as usual.

```
Tested with Charliecloud 0.37.
```

#### Other Container Runtimes

See [the Nextflow documentation for containers](https://www.nextflow.io/docs/latest/container.html).


### Configuration Files

In addition to data, `sc_pipeline.nf` expects two configuration files in JSON format: 
1. A file listing samples and associated metadata. For an example, see [Example_Input.json](/reference_files/Example_Input.json). Locations for the data files must be given as absolute paths.
2. A file setting workflow parameters. For an example, see [Example_Input_Settings.json](/reference_files/Example_Input_Settings.json). All possible parameters are listed in the `params` scope of [nextflow.config](/workflow/nextflow.config). Most parameters are optional except for `sampleInfo`, which must be the path to the sample description file described in `1`. The use of absolute paths is also recommended here, as relative paths will be interpreted relative to the execution directory.

   
### Running the Workflow

To run the workflow: 

```
nextflow run /path/to/sc_pipeline.nf [-with-report] -params-file /path/to/Input_Settings.json
```


## Input Descriptions

| Command | Description |
| --- | --- |
| sampleInfo | Full path to input file that describes the files and thier conditions (ex: reference_files/Example_Input.json) |
| clusterResolution | The leiden clustering resolution parameter (float, default: 0.4)|
| minCellsPerGene | Used to filtered on the individual samples, removes genes that do not exist in at least a number of cells (int, default: 3) |
| minGenesPerCell | Used to filtered on the individual samples, removes cells that do not have in at least a number of genes (int, default: 200)|
| mitochondrialContentMax | The maximum mitochondrial content allowed as an int, 20 means all cells with less than 20% mitochondrial content are kept (int, default: 20)|
| removeMitochondrialGenes | Removes mitochondrial genes (boolean, default: False)|
| removeRibosomalGenes | Removes ribosomal genes (boolean, default: False)|
| hvgNumTopGenes | Scanpy calculated which genes are the most highly variables and selects the top number of specified, used downstream for PCA (int, default: 5000)|
| pcaNComps |  Number of principal components selected in the principal component analysis (PCA) (int, default: 20)|
| neighborsNumber | Number of neighbors used when Scanpy calculated the nearest neighbors distance matrix and a neighborhood of graph observations (int, default: 30)|
| neighborsNpcs | Number of princpal components used when Scanpy calculated the nearest neighbors distance matrix and a neighborhood of graph observations (int, default: 20)|
| percentile | The gene expression value at this percentile is taken, all cells with an expression level greater than it are counted as having high expression for that marker gene (float, default: 90)|
| pAdj | Adjusted p-value used in the exploratory analysis (float, default: 0.05)|
| log2FC | Log2 fold change used in the exploratory analysis (float, default: 0.5)|
| heatmapNumSigGenes | The maximum number of signifant top genes used for the heatmap in the differential gene expression analysis(int, default: 50)|
| minCellsPerGroup | The number of cell required per treatment group inorder to run differential gene expression (int, default: 100)|
| nTopTermEnrichPlot | The number of top genes used in the enrichment plot (int, default: 5)|
| dotplotCutoff | The theshold for the FDR cutoff used in the enrichment dotplot (float, default: 0.25)|
| controlName | The name of the control sample group (string, default: "Control")|
| enrichTerms | An array of the datasets you would like to run gene set enrichment on, a full list of the options can be found below(string array, default: ["GO_Biological_Process_2023","GO_Molecular_Function_2023"])|
| outputFolder | Name of output folder (str, default: example_output)|
| cellType | A object with key value pairs that contain a cell type and a comma seperated string of marker genes(object key value, default: {"CD4 T cells": "CD4", "CD8 T cells": "CD8A, CD8B", "NK cells": "FCGR3A,TROBP", "B cells": "MS4A1,CD19,CD74,CD79A,IGHM", "Plasma cells": "JCHAIN,MZB1,IGHG1", "Proliferating lymphocytes: "MKI67,CD3G,FCGR3A", "Monocytes": CD14,FCGR3A,LYZ", "cDCs": "HLA-DQA1,SLC38A1", "pDCs": "BST2,MAP3K2,TRADD", "Platelets": "PF4,PPBP,ITGA2B", "Erythrocytes": HBA2,HBA1,HBB"} )|


### Full list of enrichment genesets from GSEApy
- ARCHS4_Cell-lines
- ARCHS4_IDG_Coexp
- ARCHS4_Kinases_Coexp
- ARCHS4_TFs_Coexp
- ARCHS4_Tissues
- Achilles_fitness_decrease, Achilles_fitness_increase
- Aging_Perturbations_from_GEO_down, Aging_Perturbations_from_GEO_up
- Allen_Brain_Atlas_10x_scRNA_2021, Allen_Brain_Atlas_down, Allen_Brain_Atlas_up
- Azimuth_2023
- Azimuth_Cell_Types_2021
- BioCarta_2013, BioCarta_2015, BioCarta_2016
- BioPlanet_2019
- BioPlex_2017
- CCLE_Proteomics_2020
- CORUM
- COVID-19_Related_Gene_Sets, COVID-19_Related_Gene_Sets_2021
- Cancer_Cell_Line_Encyclopedia
- CellMarker_2024
- CellMarker_Augmented_2021
- ChEA_2013, ChEA_2015, ChEA_2016, ChEA_2022
- Chromosome_Location, Chromosome_Location_hg19
- ClinVar_2019
- DGIdb_Drug_Targets_2024
- DSigDB
- Data_Acquisition_Method_Most_Popular_Genes
- DepMap_CRISPR_GeneDependency_CellLines_2023
- DepMap_WG_CRISPR_Screens_Broad_CellLines_2019
- DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019
- Descartes_Cell_Types_and_Tissue_2021
- Diabetes_Perturbations_GEO_2022
- DisGeNET
- Disease_Perturbations_from_GEO_down, Disease_Perturbations_from_GEO_up
- Disease_Signatures_from_GEO_down_2014, Disease_Signatures_from_GEO_up_2014
- DrugMatrix
- Drug_Perturbations_from_GEO_2014, Drug_Perturbations_from_GEO_down, Drug_Perturbations_from_GEO_up
- ENCODE_Histone_Modifications_2013, ENCODE_Histone_Modifications_2015
- ENCODE_TF_ChIP-seq_2014, ENCODE_TF_ChIP-seq_2015
- ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
- ESCAPE
- Elsevier_Pathway_Collection
- Enrichr_Libraries_Most_Popular_Genes
- Enrichr_Submissions_TF-Gene_Coocurrence
- Enrichr_Users_Contributed_Lists_2020
- Epigenomics_Roadmap_HM_ChIP-seq
- FANTOM6_lncRNA_KD_DEGs
- GO_Biological_Process_2013, GO_Biological_Process_2015, GO_Biological_Process_2017, GO_Biological_Process_2017b, GO_Biological_Process_2018, GO_Biological_Process_2021, GO_Biological_Process_2023
- GO_Cellular_Component_2013, GO_Cellular_Component_2015, GO_Cellular_Component_2017, GO_Cellular_Component_2017b, GO_Cellular_Component_2018, GO_Cellular_Component_2021, GO_Cellular_Component_2023
- GO_Molecular_Function_2013, GO_Molecular_Function_2015, GO_Molecular_Function_2017, GO_Molecular_Function_2017b, GO_Molecular_Function_2018, GO_Molecular_Function_2021, GO_Molecular_Function_2023
- GTEx_Aging_Signatures_2021
- GTEx_Tissue_Expression_Down, GTEx_Tissue_Expression_Up
- GTEx_Tissues_V8_2023
- GWAS_Catalog_2019, GWAS_Catalog_2023
- GeDiPNet_2023
- GeneSigDB
- Gene_Perturbations_from_GEO_down, Gene_Perturbations_from_GEO_up
- Genes_Associated_with_NIH_Grants
- Genome_Browser_PWMs
- GlyGen_Glycosylated_Proteins_2022
- HDSigDB_Human_2021
- HDSigDB_Mouse_2021
- HMDB_Metabolites
- HMS_LINCS_KinomeScan
- HomoloGene
- HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression
- HuBMAP_ASCTplusB_augmented_2022
- HumanCyc_2015, HumanCyc_2016
- Human_Gene_Atlas
- Human_Phenotype_Ontology
- IDG_Drug_Targets_2022
- InterPro_Domains_2019
- Jensen_COMPARTMENTS
- Jensen_DISEASES
- Jensen_TISSUES
- KEA_2013, KEA_2015
- KEGG_2013, KEGG_2015, KEGG_2016'
- KEGG_2019_Human, KEGG_2021_Human
- KEGG_2019_Mouse
- KOMP2_Mouse_Phenotypes_2022
- Kinase_Perturbations_from_GEO_down, Kinase_Perturbations_from_GEO_up
- L1000_Kinase_and_GPCR_Perturbations_down, L1000_Kinase_and_GPCR_Perturbations_up
- LINCS_L1000_CRISPR_KO_Consensus_Sigs
- LINCS_L1000_Chem_Pert_Consensus_Sigs
- LINCS_L1000_Chem_Pert_down, LINCS_L1000_Chem_Pert_up
- LINCS_L1000_Ligand_Perturbations_down, LINCS_L1000_Ligand_Perturbations_up
- Ligand_Perturbations_from_GEO_down, Ligand_Perturbations_from_GEO_up
- MAGMA_Drugs_and_Diseases
- MAGNET_2023
- MCF7_Perturbations_from_GEO_down, MCF7_Perturbations_from_GEO_up
- MGI_Mammalian_Phenotype_2013, MGI_Mammalian_Phenotype_2017, MGI_Mammalian_Phenotype_Level_3, MGI_Mammalian_Phenotype_Level_4, MGI_Mammalian_Phenotype_Level_4_2019, MGI_Mammalian_Phenotype_Level_4_2021, MGI_Mammalian_Phenotype_Level_4_2024
- MSigDB_Computational
- MSigDB_Hallmark_2020
- MSigDB_Oncogenic_Signatures
- Metabolomics_Workbench_Metabolites_2022
- Microbe_Perturbations_from_GEO_down, Microbe_Perturbations_from_GEO_up
- MoTrPAC_2023
- Mouse_Gene_Atlas
- NCI-60_Cancer_Cell_Lines
- NCI-Nature_2016
- NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions, NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions, NIH_Funded_PIs_2017_Human_AutoRIF, NIH_Funded_PIs_2017_Human_GeneRIF
- NURSA_Human_Endogenous_Complexome
- OMIM_Disease
- OMIM_Expanded
- Old_CMAP_down, Old_CMAP_up
- Orphanet_Augmented_2021
- PFOCR_Pathways, PFOCR_Pathways_2023
- PPI_Hub_Proteins
- PanglaoDB_Augmented_2021
- Panther_2015, Panther_2016
- Pfam_Domains_2019
- Pfam_InterPro_Domains
- PheWeb_2019
- PhenGenI_Association_2021
- Phosphatase_Substrates_from_DEPOD
- ProteomicsDB_2020
- Proteomics_Drug_Atlas_2023
- RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO
- RNAseq_Automatic_GEO_Signatures_Human_Down, RNAseq_Automatic_GEO_Signatures_Human_Up
- RNAseq_Automatic_GEO_Signatures_Mouse_Down, RNAseq_Automatic_GEO_Signatures_Mouse_Up
- Rare_Diseases_AutoRIF_ARCHS4_Predictions, Rare_Diseases_AutoRIF_Gene_Lists, Rare_Diseases_GeneRIF_ARCHS4_Prediction, Rare_Diseases_GeneRIF_Gene_Lists
- Reactome_2013, Reactome_2015, Reactome_2016, Reactome_2022
- Rummagene_kinases
- Rummagene_signatures
- Rummagene_transcription_factors
- SILAC_Phosphoproteomics
- SubCell_BarCode
- SynGO_2022, SynGO_2024
- SysMyo_Muscle_Gene_Sets
- TF-LOF_Expression_from_GEO
- TF_Perturbations_Followed_by_Expression
- TG_GATES_2020
- TRANSFAC_and_JASPAR_PWMs
- TRRUST_Transcription_Factors_2019
- Table_Mining_of_CRISPR_Studies
- Tabula_Muris
- Tabula_Sapiens
- TargetScan_microRNA, TargetScan_microRNA_2017
- The_Kinase_Library_2023
- Tissue_Protein_Expression_from_Human_Proteome_Map
- Tissue_Protein_Expression_from_ProteomicsDB
- Transcription_Factor_PPIs
- UK_Biobank_GWAS_v1
- Virus-Host_PPI_P-HIPSTer_2020
- VirusMINT
- Virus_Perturbations_from_GEO_down, Virus_Perturbations_from_GEO_up
- WikiPathway_2021_Human, WikiPathway_2023_Human
- WikiPathways_2013, WikiPathways_2015, WikiPathways_2016
- WikiPathways_2019_Human
- WikiPathways_2019_Mouse
- dbGaP
- huMAP
- lncHUB_lncRNA_Co-Expression
- miRTarBase_2017
