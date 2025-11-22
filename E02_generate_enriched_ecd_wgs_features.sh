while getopts "i:o:f:n:s:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    f )
      path_to_fa=$OPTARG
      ;;
    n )
      path_to_nucleosome_ref=$OPTARG
      ;;
    s )
      feature_srcdir=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-r] run_name"
      exit 1
      ;;
  esac
done

sampleid=$(echo ${inputbam} | xargs -n 1 basename);
sampleid=$(echo ${sampleid} | cut -d '.' -f 1)
mkdir -p ${outputdir}/${sampleid};

echo -e "WORKING ON BAM FILE " $sampleid " with bed file " $run_name;

##### Step 1: convert BAM file to FRAG file.
if [ ! -f ${outputdir}/${sampleid}.frag.tsv ]; then
    bash ${feature_srcdir}/01_convert_BAM_to_FRAGS.sh \
        -i ${inputbam} \
        -o ${outputdir} \
        -n 20 \
        -f ${path_to_fa} \
        -q 50 \
        -w 150 \
        -e 151 \
        -r 350 \
        -c "true"
    echo -e "finished converting BAM to FRAG for sample " ${sampleid};
else
    echo -e "BAM to FRAG for sample " ${sampleid} " with strategy " ${strategy} " is done!\n";
    echo -e "output files are saved to " ${outputdir}"\n";
fi


##### Step 2: generate features from FRAG file.
if [ ! -f ${outputdir}/${sampleid}.frag_output.tsv ]; then
    bash ${feature_srcdir}/02_calculate_features_from_frag.sh \
        -i ${outputdir}/${sampleid}/${sampleid}.frag.tsv \
        -o ${outputdir} \
        -f ${path_to_fa} \
        -r ${path_to_nucleosome_ref} \
        -c "true"
else
    echo -e "Feature calculation for sample " ${sampleid} " with strategy " ${strategy} " is done!\n";
    echo -e "output files are saved to " ${outputdir}"\n";
fi

##### Step 3: generate final feature files.
echo -e "Generating final feature files ..."
# source /home/hieunguyen/miniconda3/bin/activate && conda activate pytorch;
echo ${outputdir}/${sampleid}/${sampleid}.frag_output.tsv
echo  ${outputdir}/${sampleid}/fragmentomics_features
python ${feature_srcdir}/04_generate_features_fragmentomics_and_image_features.py \
    --input ${outputdir}/${sampleid}/${sampleid}.frag_output.tsv \
    --output ${outputdir}/${sampleid}/fragmentomics_features;

# finished!!!!!