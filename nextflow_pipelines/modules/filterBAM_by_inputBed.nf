process filterBAM_by_inputBed {
    // filter a BAM file by a input BED file
    tag "$sample_id"
    cache "deep";
    
    input:
        tuple val(sample_id), file(inputbam)
        file(inputbed)
        file(inputSH)
    output:
        tuple val(sample_id), file("*${sample_id}*.filtered.bam"), emit: filtered_bam

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bash ${inputSH} -i ${inputbam} -b ${inputbed} -o . 
    """
}