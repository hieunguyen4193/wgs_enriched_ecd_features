while getopts "i:b:o:r:" opt; do
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
    r )
      run_name=$OPTARG
      ;;
    
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-r] run_name"
      exit 1
      ;;
  esac
done

# we use run_name to mark the name of bed file reigons!
mkdir -p ${outputdir};
mkdir -p ${outputdir}/${run_name};

echo -e "output files will be saved at ${outputdir}\n"

sampleid=$(echo ${inputbam} | xargs -n 1 basename);
echo -e "working on file "  $sampleid;
samtools view -b -h -L ${inputbed} ${inputbam} > ${outputdir}/${run_name}/${sampleid%.bam*}.filtered.bam;
samtools index ${outputdir}/${run_name}/${sampleid%.bam*}.filtered.bam;