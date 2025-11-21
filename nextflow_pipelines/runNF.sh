
nextflow run main.nf \
    --SAMPLE_SHEET $SAMPLE_SHEET \
    --INPUT_BED $INPUT_BED \
    --OUTDIR $OUTDIR \
    --E01_SH $E01_SH \
    --E02_SH $E02_SH \
    --E03_SH $E03_SH \
    --RUN_NAME $RUN_NAME \
    --PATH_TO_FA $PATH_TO_FA \
    --PANDEPTH $PANDEPTH \
    -resume -c default.config \
    -w ${WORKDIR} \
    -with-report "${OUTDIR}/report.html" \
    -with-timeline "${OUTDIR}/timeline.html" \
    -with-dag "${OUTDIR}/dag.svg";
