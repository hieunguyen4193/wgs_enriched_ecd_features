workflow pipeline_init {
    take:
    //  path to the input samplesheet .csv file, the sampleshet file should contain the columns SampleID, BAM file
    //  other format will not be accepted. 
    // Please modify your sampleshete file accordingly.
        input_SampleSheet 
    main:
        Channel
        .fromPath(input_SampleSheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.SampleID, file(row.BAM_path))}
        .set{input_bam_ch}
    
    emit:
    samplesheet = input_bam_ch // emit to the samplesheet channel, use as input for other downstream processes
}