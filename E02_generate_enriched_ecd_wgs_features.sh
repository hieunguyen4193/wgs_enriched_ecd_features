while getopts "i:o:r:f:" opt; do
  case ${opt} in
    i )
      inputbam=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    r )
      run_name=$OPTARG
      ;;
    f )
      path_to_fa=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-i] inputbam [-o] outputdir [-r] run_name"
      exit 1
      ;;
  esac
done

sampleid=$(echo ${inputbam} | xargs -n 1 basename);
mkdir -p ${outputdir}/${run_name}/${sampleid};

echo -e "WORKING ON BAM FILE " $sampleid " with bed file " $run_name;

##### Step 1: convert BAM file to FRAG file.
if [ ! -f ${outputdir}/${run_name}/${sampleid}.frag.tsv ]; then
    bash ${feature_srcdir}/01_convert_BAM_to_FRAGS.sh \
        -i ${inputbam} \
        -o ${outputdir}/${run_name} \
        -n 20 \
        -f ${path_to_fa} \
        -q 50 \
        -w 150 \
        -e 151 \
        -r 350 \
        -c "true"
else
    echo -e "BAM to FRAG for sample " ${sampleid} " with strategy " ${strategy} " is done!\n";
    echo -e "output files are saved to " ${outputdir}/${run_name}"\n";
fi


##### Step 2: generate features from FRAG file.
if [ ! -f ${outputdir}/${run_name}/${sampleid}.frag_output.tsv ]; then
    bash ${feature_srcdir}/02_calculate_features_from_frag.sh \
        -i ${outputdir}/${run_name}/${sampleid}/${sampleid}.frag.tsv \
        -o ${outputdir}/${run_name} \
        -f ${path_to_fa} \
        -r ${path_to_nucleosome_ref} \
        -c "true"
else
    echo -e "Feature calculation for sample " ${sampleid} " with strategy " ${strategy} " is done!\n";
    echo -e "output files are saved to " ${outputdir}/${run_name}"\n";
fi

##### Step 3: generate final feature files.
echo -e "Generating final feature files ..."
# source /home/hieunguyen/miniconda3/bin/activate && conda activate pytorch;
echo ${outputdir}/${run_name}/${sampleid}/${sampleid}.frag_output.tsv
echo  ${outputdir}/${run_name}/${sampleid}/fragmentomics_features
python ${feature_srcdir}/04_generate_features_fragmentomics_and_image_features.py \
    --input ${outputdir}/${run_name}/${sampleid}/${sampleid}.frag_output.tsv \
    --output ${outputdir}/${run_name}/${sampleid}/fragmentomics_features;

# finished!!!!!