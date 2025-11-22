process generate_enriched_ecd_wgs_features {
    // filter a BAM file by a input BED file
    tag "$sample_id"
    cache "deep";
    
    input:
        tuple val(sample_id), file(inputbam)
        file(inputSH)
        val(run_name)
        val(path_to_fa)
        val(feature_srcdir)
    output:
        tuple val(sample_id), file("${sample_id}.filtered.bam")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bash ${inputSH} -i ${inputbam} -o . -r ${run_name} -f ${path_to_fa} -s {feature_srcdir}
    """
}