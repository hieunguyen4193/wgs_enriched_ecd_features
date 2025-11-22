while getopts "i:b:o:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    b )
      inputbed=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir"
      exit 1
      ;;
  esac
done

mkdir -p ${outputdir};

echo -e "output files will be saved at ${outputdir}\n"

sampleid=$(echo ${inputbam} | xargs -n 1 basename);
echo -e "working on file "  $sampleid;
samtools view -b -h -L ${inputbed} ${inputbam} > ${outputdir}/${sampleid%.bam*}.filtered.bam;
samtools index ${outputdir}/${sampleid%.bam*}.filtered.bam;