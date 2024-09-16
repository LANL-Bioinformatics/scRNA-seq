#!/usr/bin/env nextflow

import groovy.json.JsonSlurper

//takes qc-filtered h5ad files. 
//Produces a series of graphs and integrated.h5ad file with all samples combined, normalized, and batch corrected
process integrated_dataset_qc {
    publishDir(
        path: "${params.outputFolder}",
        mode: 'copy'
    )

    input:
    path annData

    output:
    path "*.png"
    path "integrated.h5ad"

    script:
    //setup optional arguments
    def hvgTop = params.hvgNumTopGenes != null ? "--HVG_num_top_genes $params.hvgNumTopGenes" : ""
    def pcaComps = params.pcaNcomps != null ? "--PCA_N_Comps $params.pcaNcomps" : ""
    def nNeighbors = params.neighborsNumber != null ? "--Neighbors_Number $params.neighborsNumber" : ""
    def nNeighborPCs = params.neighborsNpcs != null ? "--Neighbors_N_PCS $params.neighborsNpcs" : ""

    //invoke script
    """
    ${workflow.projectDir}/../scripts/QC_intergrated_dataset.py \
    $hvgTop \
    $pcaComps \
    $nNeighbors \
    $nNeighborPCs \
    --in_files $annData \
    --output_folder \$PWD
    """
}

//takes qc-filtered h5ad files.
//Produces a series of plots comparing the samples inside each treatment/outcome group.
process plots_across_groups_qc {
    publishDir(
        path: "${params.outputFolder}",
        mode: 'copy'
    )
    input:
    path annData

    output:
    path "*.png"

    script:
    //invoke script
    """
    ${workflow.projectDir}/../scripts/QC_plots_across_groups.py \
    --in_files $annData \
    --output_folder \$PWD \
    """
}

//takes in one sample (h5 or mtx/tsv files) and associated metadata.
//Runs sample-specific QC and outputs h5ad files that preserve annotation
process per_sample_qc {
    publishDir(
        path: "${params.outputFolder}",
        mode: 'copy'
    )
    input:
    tuple val(md), path(normalFile), path(rawFile), path(ribosomeList)

    output:
    path "*.png" //plots
    path "*.h5ad", emit: annData

    script:
    //null checks and defining arguments
    def sampleName = md["sampleName"] != null ? "--sample_name \"${md["sampleName"]}\"" : ""
    def condition = md["condition"] != null ? "--condition \"${md["condition"]}\"" : ""
    def batch = md["batch"] != null ? "--batch \"${md["batch"]}\"" : ""

    def clusterRes = params.clusterResolution != null ? "--cluster_resolution ${params.clusterResolution}" : ""
    def minCells = params.minCellsPerGene != null ? "--min_cells_per_gene ${params.minCellsPerGene}" : ""
    def minGenes = params.minGenesPerCell != null ? "--min_genes_per_cell ${params.minGenesPerCell}" : ""
    def mitMax = params.mitochondrialContentMax != null ? "--mitochondrial_content_max ${params.mitochondrialContentMax}" : ""


    //For optional inputs (since 1 or the other type of files can be provided):
    //Only generate the CL argument if a file of this type was provided, and generate the correct argument
    //for the files provided. A list of files (eg, mtx_tsv files) will expand in interpolation 

    def regFileArg = ""
    def rawFileArg = ""
    
    if(normalFile.name != "${workflow.projectDir}/nf_assets/NO_FILE") {
        if(md["regFormat"] == "h5") {
            regFileArg = "--h5_file_name ${normalFile}"
        } else {
            regFileArg = "--mtx_tsv_file_name ${normalFile}"
        }
    }

    if(rawFile.name != "${workflow.projectDir}/nf_assets/NO_FILE2") {
        if(md["rawFormat"] == "h5") {
            rawFileArg = "--h5_raw_file_name \"${rawFile}\""
        } else {
            rawFileArg = "--mtx_tsv_raw_file_name \"${rawFile}\""
            
        }
    }

    
    //invoke script with arguments
    """
    ${workflow.projectDir}/../scripts/QC_per_sample.py \
    $sampleName \
    $condition \
    $batch \
    $regFileArg \
    $rawFileArg \
    $clusterRes \
    $minCells \
    $minGenes \
    $mitMax \
    --ribosomal_genelist $ribosomeList \
    --output_folder \$PWD
    """
}

workflow {
    "mkdir ${workflow.projectDir}/nf_assets".execute().text
    "touch ${workflow.projectDir}/nf_assets/NO_FILE".execute().text
    "touch ${workflow.projectDir}/nf_assets/NO_FILE2".execute().text

    //parse input
    def jsonSlurper = new JsonSlurper()
    def info = jsonSlurper.parse(new File(params.sampleInfo))
    
    //We're separating out the file paths, since Nextflow needs to be able to stage those in its work directories
    def regFileList = []
    def rawFileList = []
    def metadataList = []
    for (sample in info.Samples) {
        rawType = ""
        regType = ""
        //separate out raw file paths
        if(sample["Raw_File_Name"] != null) {
            rawFileList.add(new ArrayList(sample["Raw_File_Name"].values()))
            if(rawFileList[-1][0].endsWith(".h5")) { //adding some information for downstream processes
                rawType = "h5"
            } else {
                rawType = "mtx_tsv"
            }
        } else {
            rawFileList.add(["${workflow.projectDir}/nf_assets/NO_FILE"])
        }

        //separate out normal file paths
        if(sample["File_Name"] != null) {
            regFileList.add(new ArrayList(sample["File_Name"].values()))
            if(regFileList[-1][0].endsWith(".h5")) {
                regType = "h5"
            } else {
                regType = "mtx_tsv"
            }
        } else {
            regFileList.add(["${workflow.projectDir}/nf_assets/NO_FILE2"])
        }

        //separate out associated metadata
        metadataList.add(
            [sampleName: sample["Sample_Name"],  //elvis operate these for nullchecks?
            condition: sample["Condition"],
            batch: sample["Batch"],
            rawFormat: rawType,
            regFormat: regType])
    }
    //current state: three lists with matching entries for the metadata, regular file path(s), and raw file path(s)
    def tupleList = []
    for(int i=0;i<metadataList.size();i++) {
        tupleList.add([metadataList[i],regFileList[i],rawFileList[i]])
    }
    //new state: one list of tuples with three elements: metadata, regular file path(s), and raw file path(s)

    //combine that list of inputs with a path to the file with a list of ribosomal genes
    //necessary to specify that file in this workflow so that it can be staged into the working directory
    input_ch = channel.fromList(tupleList).combine(channel.fromPath(params.ribosomalGeneList, checkIfExists:true))

    //invoke per-sample qc script
    per_sample_qc(input_ch)
    qc_filt_ch = per_sample_qc.out.annData.collect()

    //run both other qc scripts (plots and integrated QC) in parallel off of the output from the per-sample QC
    plots_across_groups_qc(qc_filt_ch)
    integrated_dataset_qc(qc_filt_ch)


}