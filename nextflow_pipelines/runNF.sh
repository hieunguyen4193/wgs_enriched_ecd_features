INPUT_BED="/mnt/nvme/DATA_HIEUNGUYEN/resource/assets/TSS_UTR5p_regions/04_output/Normal/TSS_beds/Breast/v0.1/down_genes/promoter_regions_up_1000_down_1000.bed"
RUN_NAME="down_genes_v0.1_Breast_TSS_beds_Normal";

SAMPLE_SHEET="SampleSheet_BAM_highdepth.csv";
MAIN_OUTDIR="/mnt/nvme/DATA_HIEUNGUYEN/outdir/TOO_enriched_features";
OUTDIR="${MAIN_OUTDIR}/${RUN_NAME}";
mkdir -p $OUTDIR;

E01_SH="../E01_generate_filter_BAM_files.sh";
E02_SH="../E02_generate_enriched_ecd_wgs_features.sh";
E03_SH="../E03_generate_TSS_count_features_at_enriched_regions.sh"
PATH_TO_FA="/mnt/nvme/DATA_HIEUNGUYEN/resource/hg19.fa"

# for feature_srcdir download the ecd_wgs_feature pipeline from github - release version
# https://github.com/gsecddatalab/ecd_wgs_feature_pipeline/releases/tag/ecd_wgs_pipeline_v1.0

feature_srcdir="./ecd_wgs_feature_pipeline-ecd_wgs_pipeline_v1.0"

WORKDIR="/mnt/nvme/DATA_HIEUNGUYEN/outdir/TOO_enriched_features/workdir"
mkdir -p $WORKDIR;
PANDEPTH="/mnt/nvme/DATA_HIEUNGUYEN/resource/PanDepth/bin/pandepth";

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
    --feature_srcdir $feature_srcdir \
    -resume -c default.config \
    -w ${WORKDIR} \
    -with-report "${OUTDIR}/report.html" \
    -with-timeline "${OUTDIR}/timeline.html" \
    -with-dag "${OUTDIR}/dag.svg";
