process generate_enriched_TSS_count_features {
    // filter a BAM file by a input BED file
    tag "$sample_id"
    cache "deep";
    
    input:
        tuple val(sample_id), file(inputbam)
        file(inputSH)
        file(inputbed)
        val(run_name)
        val(pandepth)
    output:
        tuple val(sample_id), file("${sample_id}.filtered.bam")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bash ${inputSH} -i ${inputbam} -b ${inputbed} -o . -p ${pandepth} -r ${run_name}
    """
}