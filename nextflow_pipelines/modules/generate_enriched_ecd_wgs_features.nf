process generate_enriched_ecd_wgs_features {
    // filter a BAM file by a input BED file
    tag "$sample_id"
    cache "deep";
    
    input:
        tuple val(sample_id), file(inputbam)
        file(inputSH)
        val(path_to_fa)
        file(feature_srcdir)
        file(path_to_nucleosome_ref)
    output:
        tuple val(sample_id), file("*${sample_id}*")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bash ${inputSH} -i ${inputbam} -o . -f ${path_to_fa} -s ${feature_srcdir} -n ${path_to_nucleosome_ref}
    """
}