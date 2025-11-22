
include { filterBAM_by_inputBed as filterBAM_by_inputBed} from "../modules/filterBAM_by_inputBed.nf"
include { generate_enriched_ecd_wgs_features } from "../modules/generate_enriched_ecd_wgs_features.nf"
include { generate_enriched_TSS_count_features } from "../modules/generate_enriched_TSS_count_features.nf"

workflow run_enriched_features_and_TSS_count_features {
    take:
        inputbam
        inputbed
        E01_sh
        E02_sh
        E03_sh
        run_name
        path_to_fa
        pandepth
        feature_srcdir
    main:
        filtered_bam = filterBAM_by_inputBed(inputbam, inputbed, E01_sh, run_name)
        generate_enriched_ecd_wgs_features(filtered_bam, E02_sh, run_name, path_to_fa, feature_srcdir)
        generate_enriched_TSS_count_features(filtered_bam, E03_sh, inputbed, run_name, pandepth)
}   