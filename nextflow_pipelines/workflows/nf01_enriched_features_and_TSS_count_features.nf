include { run_enriched_features_and_TSS_count_features } from "../subworkflows/run_enriched_features_and_TSS_count_features.nf"
include { pipeline_init }  from "../subworkflows/pipeline_init.nf"

workflow nf01_enriched_features_and_TSS_count_features {
    take:
        input_SampleSheet // path to the input samplesheet .csv file, the sampleshet file should contain the columns SampleID, FASTQ1, and FASTQ2
        inputbed
        E01_sh
        E02_sh
        E03_sh
        path_to_fa
        pandepth
        feature_srcdir
        path_to_nucleosome_ref
    main:
        pipeline_init(input_SampleSheet)   
        run_enriched_features_and_TSS_count_features(
            pipeline_init.out.samplesheet,
            inputbed,
            E01_sh,
            E02_sh,
            E03_sh,
            path_to_fa,
            pandepth,
            feature_srcdir,
            path_to_nucleosome_ref
        )
}