# this script is tailored for running in the GS-HPC only. 
# please modify accordingly if you want to run in other machine
#####

mkdir -p /mnt/DATASM14/DATA_HIEUNGUYEN/outdir/enriched_features/tmp_for_sort
export PATH=/home/hieunguyen/samtools/bin:$PATH;
export TMPDIR=/mnt/DATASM14/DATA_HIEUNGUYEN/outdir/enriched_features/tmp_for_sort
maindir="/mnt/DATASM14/DATA_HIEUNGUYEN/src/ecd_wgs_enriched_and_kmer_features";

bash ${maindir}/enrichedFeature_pipeline/atacseq_bed/filter_BAM_files_and_generate_features.sh;

# bash ${maindir}/enrichedFeature_pipeline/bin600/filter_BAM_files_and_generate_features.sh;

# bash ${maindir}/enrichedFeature_pipeline/custom_genes_TSS/filter_BAM_files_and_generate_features.sh;

# bash ${maindir}/enrichedFeature_pipeline/methylation_regions/filter_BAM_files_and_generate_features.sh;

# bash ${maindir}/enrichedFeature_pipeline/panel_7_8_vs_nucleosomeMap/filter_BAM_files_and_generate_features.sh;

# bash ${maindir}/enrichedFeature_pipeline/panel_7_8_vs_TSS/filter_BAM_files_and_generate_features.sh;

# bash ${maindir}/enrichedFeature_pipeline/TCGA_atac_seq_beds/filter_BAM_files_and_generate_features.sh;