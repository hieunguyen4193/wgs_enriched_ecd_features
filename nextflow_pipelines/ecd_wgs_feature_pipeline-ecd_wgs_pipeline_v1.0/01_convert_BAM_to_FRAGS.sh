#####----------------------------------------------------------------------#####
##### INTRODUCTION
#####----------------------------------------------------------------------#####
# This script pre-process an input BAM file to a 
# fragment-wise data features, which can be use to calculate
# several fragmentomics features. 
# export PATH=/Users/hieunguyen/samtools/bin:$PATH
# bash 01_BAM_to_FRAGS.sh -i /Users/hieunguyen/data/bam_files/WGShg19.bam  -o ./output/ -n 10

#####----------------------------------------------------------------------#####
##### input args
#####----------------------------------------------------------------------#####
while getopts "i:o:n:f:q:w:e:r:c:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    n )
      samtools_num_threads=$OPTARG
      ;;
    f )
      ref=$OPTARG
      ;;
    q )
      short_lower=$OPTARG
      ;;
    w )
      short_upper=$OPTARG
      ;;
    e )
      long_lower=$OPTARG
      ;;
    r )
      long_upper=$OPTARG
      ;;
    c )
      cleanup=$OPTARG
      ;;
       
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads [-f] reference genome"
      exit 1
      ;;
  esac
done


echo -e "input bam file: " ${inputbam}
# Check if the input BAM file exists
if [ ! -f "${inputbam}" ]; then
    echo "Input BAM file does not exist: ${inputbam}"
    exit 1
fi

# Generate the output folder for the input sample
sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=$(echo ${sampleid} | cut -d '.' -f 1)
outputdir=${outputdir}/${sampleid}

mkdir -p ${outputdir}

# picard="./picard.jar";
picard="/home/hieunguyen/src/wgs_enriched_ecd_features/picard.jar"
# get picard from wget https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

# preprocessing the input BAM file, check if the file exists, if yes, run the script, if yes, skip
if [ ! -f "${outputdir}/${sampleid}.sortedN.markdup.bam" ]; then
    echo -e "Pre-processing the input BAM file " ${sampleid} " ..."
    ## sort base on read name and remove unpaired and unmapped reads from BAM files, sort BAM files by read name,
    samtools view -h -f 3 ${inputbam} | samtools sort -n -@ ${samtools_num_threads} -o ${outputdir}/${sampleid}.sortedN.bam;

    ## mark duplicates BAM file with picard
    # java -Xms512m -Xmx4g -jar ./picard.jar MarkDuplicates \
    java -Xms512m -Xmx4g -jar ${picard} MarkDuplicates \
        -I ${outputdir}/${sampleid}.sortedN.bam \
        -O ${outputdir}/${sampleid}.sortedN.markdup.bam \
        -M ${outputdir}/${sampleid}.marked_dup_metrics.txt
else
    echo -e "The input BAM file " ${sampleid} " has been pre-processed and marked duplicates. Skip this step ..."
fi

# remove flen == 0 reads, keep only reads mapped to chr1-22, remove all others
#  and generate information on mate reads, fragment wise level data table
# finall, sort file based on read name and read position coordinates
        
if [ ! -f "${outputdir}/${sampleid}.frag.tsv" ]; then
    echo -e "Converting the input BAM file to FRAG files ..."
    samtools view ${outputdir}/${sampleid}.sortedN.markdup.bam | \
    awk -v OFS='\t' '{if ($9 != 0 && $3 ~ /^chr([1-9]|1[0-9]|2[0-2])$/){print $1 "\t" $3 "\t" $4 "\t" $8 "\t" $9 "\t" $5}}' | \
    awk -v OFS='\t' '{if ($5 > 0){$7=$3+$5; $8=$1"_"$2"_"$3; print $0} else {$7=$4-$5; $3=$4; $8=$1"_"$2"_"$3; print $0} }' | \
    sort -k8,8 | \
    awk '{ print $2 "\t" $3 "\t" $7 "\t" $5 "\t" $8 "\t" $6}' > ${outputdir}/${sampleid}.frag.tsv
else 
    echo -e "The input BAM file " ${sampleid} " has been converted to FRAG files. Skip this step ..."
fi

# split BAM file to short and long fragments BAM
if [ ! -f "${outputdir}/${sampleid}.finished_splitBAM.txt" ]; then
    ## split BAM file to .short.bam and .long.bam files
    echo -e "Splitting the input BAM file into short and long fragments ..."

    samtools sort ${outputdir}/${sampleid}.sortedN.markdup.bam -@ ${samtools_num_threads} -o ${outputdir}/${sampleid}.markdup.sorted.bam
    bash split_BAM_to_short_and_long.sh \
        -i ${outputdir}/${sampleid}.markdup.sorted.bam \
        -o ${outputdir}/short_long_BAM \
        -n ${samtools_num_threads} \
        -q ${short_lower} \
        -w ${short_upper} \
        -e ${long_lower} \
        -r ${long_upper}
    touch ${outputdir}/${sampleid}.finished_splitBAM.txt
else
    echo -e "The input BAM file " ${sampleid} " has been split into short and long fragments. Skip this step ..."
fi

# clean-up
if [ "${cleanup}" == "true" ]; then
    echo "Cleaning up intermediate files..."
    rm -rf ${outputdir}/${sampleid}.markdup.sorted.bam
    rm -rf ${outputdir}/${sampleid}.markdup.sorted.bam.bai
fi