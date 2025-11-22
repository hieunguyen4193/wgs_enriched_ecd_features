export PATH=/home/hieunguyen/bedtools2/bin:$PATH

inputbam="/media/hieunguyen/GSHD_HN01/raw_data/bam_files/WGShg19.bam";
path_to_fa="/media/hieunguyen/GSHD_HN01/storage/resources/hg19.fa";
path_to_nucleosome_ref="/media/hieunguyen/GSHD_HN01/storage/resources/rpr_map_EXP0779.sorted.bed";
outputdir="./output";

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=$(echo ${sampleid} | cut -d '.' -f 1)

bash 01_convert_BAM_to_FRAGS.sh -i ${inputbam} -o ${outputdir} -n 20 -f ${path_to_fa} -q 50 -w 150 -e 151 -r 350 -c "true"

bash 02_calculate_features_from_frag.sh -i ${outputdir}/${sampleid}/${sampleid}.frag.tsv -o ${outputdir} -f ${path_to_fa} -r ${path_to_nucleosome_ref} -c "true"

# this Rscript needs to be run in DOCKER: tronghieunguyen/ecd_features
# Rscript 03_generate_binWise_features.R \
#     -s ${outputdir}/${sampleid}/short_long_BAM/${sampleid}.markdup.sorted_50_150.short.bam \
#     -l ${outputdir}/${sampleid}/short_long_BAM/${sampleid}.markdup.sorted_50_150.long.bam \
#     -f ${outputdir}/${sampleid}/${sampleid}.markdup.sorted.bam \
#     -o ${outputdir}/${sampleid}/binwise_features;

echo -e "Generating final feature files ..."
source /home/hieunguyen/miniconda3/bin/activate && conda activate pytorch;
echo ${outputdir}/${sampleid}/${sampleid}.frag_output.tsv
echo  ${outputdir}/${sampleid}/fragmentomics_features
python 04_generate_features_fragmentomics_and_image_features.py \
    --input ${outputdir}/${sampleid}/${sampleid}.frag_output.tsv \
    --output ${outputdir}/${sampleid}/fragmentomics_features;

