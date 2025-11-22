while getopts "i:o:n:t:q:w:e:r:" opt; do
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
    
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-n] samtools_num_threads -[q] short_lower -[w] short_upper -[e] long_lower -[r] long_upper"
      exit 1
      ;;
  esac
done
mkdir -p ${outputdir};
echo -e "input bam file: " ${inputbam}
# Check if the input BAM file exists
if [ ! -f "${inputbam}" ]; then
    echo "Input BAM file does not exist: ${inputbam}"
    exit 1
fi

sampleid=$(echo ${inputbam} | xargs -n 1 basename)
sampleid=${sampleid%.bam*}

echo -e "splitting BAM file into short and long BAM files ..."

# split BAM file: short_lower <= insert size <= short_upper
samtools view -h ${inputbam} | \
    awk 'substr($0,1,1)=="@" || ($9 >= '${short_lower}' && $9 <= '${short_upper}') || ($9 <= -'${short_lower}' && $9 >= -'${short_upper}')' | \
    samtools view -b | samtools sort -@ ${samtools_num_threads} > ${outputdir}/${sampleid}_${short_lower}_${short_upper}.short.bam;

# split BAM file: long_lower <= insert size <= long_upper
samtools view -h ${inputbam} | \
    awk 'substr($0,1,1)=="@" || ($9 >= '${long_lower}' && $9 <= '${long_upper}') || ($9 <= -'${long_lower}' && $9 >= -'${long_upper}')' | \
    samtools view -b | samtools sort -@ ${samtools_num_threads} > ${outputdir}/${sampleid}_${long_lower}_${long_upper}.long.bam;

samtools index ${outputdir}/${sampleid}_${short_lower}_${short_upper}.short.bam;
samtools index ${outputdir}/${sampleid}_${long_lower}_${long_upper}.long.bam;

touch ${outputdir}/${sampleid}.finished_splitBAM.txt