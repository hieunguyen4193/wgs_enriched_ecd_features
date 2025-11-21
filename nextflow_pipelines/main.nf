#!/usr/bin/env nextflow
/*
NEXTFLOW DSL 2 pipeline for single cell data analysis
This pipeline is designed to process single-cell RNA sequencing data using the 10X Genomics Cell Ranger pipeline.

tronghieunguyen@pm.me
*/

// enable nextflow DSL2

nextflow.enable.dsl = 2

include { nf01_enriched_features_and_TSS_count_features } from "./workflows/nf01_enriched_features_and_TSS_count_features.nf"

workflow {
    main:
    //  run the main pipieline. 
    nf01_enriched_features_and_TSS_count_features(
        file(params.SAMPLE_SHEET), 
        file(params.INPUT_BED),
        file(params.E01_SH),
        file(params.E02_SH),
        file(params.E03_SH),
        params.RUN_NAME,
        params.PATH_TO_FA,
        params.PANDEPTH
    )
}