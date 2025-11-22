while getopts "i:o:f:r:c:" opt; do
  case ${opt} in
    i )
      inputfrag=$OPTARG
      ;;
    o )
      outputdir=$OPTARG
      ;;
    f )
      path_to_fa=$OPTARG
      ;;
    r )
      nucleosome_ref=$OPTARG
      ;;  
    c )
      cleanup=$OPTARG
      ;;  
    \? )
      echo "Usage: cmd [-i] input .frag file [-o] outputdir [-f] path_to_fa [-r] nucleosome_ref [-n] NDR TOO ref [-b] NDR binary ref [-c] cleanup"
      exit 1
      ;;
  esac
done

sampleid=$(echo ${inputfrag} | xargs -n 1 basename)
sampleid=${sampleid%.frag*}
outputdir=${outputdir}/${sampleid}
mkdir -p ${outputdir};

echo -e "-------------------------------------------------------------------------"
echo -e "Generating features from fragment file " ${inputfrag} " ..."
echo -e "-------------------------------------------------------------------------"
#####----------------------------------------------------------------------#####
##### 4bp END MOTIF
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.finished_4bpEM.txt" ]; then
    echo -e "getting 4bp end motif"
    cat ${inputfrag} | \
      awk '{start=$2 - 1; end= $2 - 1 + 4; name= $5; strand = "+"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.forward_endcoord4bp.bed;

    cat ${inputfrag} | \
      awk '{start=$3 - 1 - 4; end= $3 - 1; name= $5; strand = "-"; print $1 "\t" start "\t" end "\t" name "\t" "1" "\t" strand}' \
    > ${outputdir}/${sampleid}.reverse_endcoord4bp.bed;

    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.forward_endcoord4bp.bed | \
      awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${sampleid}.forward_endmotif4bp.txt
    bedtools getfasta -s -name -tab -fi ${path_to_fa} -bed ${outputdir}/${sampleid}.reverse_endcoord4bp.bed | \
    awk -v OFS='\t' '{split($0, a, "::"); $1=a[1]; print $0}'  > ${outputdir}/${sampleid}.reverse_endmotif4bp.txt

    rm -rf ${outputdir}/${sampleid}.forward_endcoord4bp.bed
    rm -rf ${outputdir}/${sampleid}.reverse_endcoord4bp.bed

    sort -k1,1 ${outputdir}/${sampleid}.forward_endmotif4bp.txt > ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt
    sort -k1,1 ${outputdir}/${sampleid}.reverse_endmotif4bp.txt > ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt

    touch ${outputdir}/${sampleid}.finished_4bpEM.txt
fi

count_4bpEM_forward=$(cat ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt | wc -l)
count_4bpEM_reverse=$(cat ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt | wc -l)

#####----------------------------------------------------------------------#####
##### NUCLEOSOME FOOTPRINT
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.finished_Nucleosome.txt" ]; then
  echo -e "generating nucleosome features ..."
  cat ${inputfrag} | cut -f1,2,4,5 | \
    awk -v OFS='\t' '{$5=$2 + 1; print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.forward_Nucleosome.bed
  cat ${inputfrag} | cut -f1,3,4,5 | \
    awk -v OFS='\t' '{$5=$2 + 1; print $1 "\t" $2 "\t" $5 "\t" $4 "\t" $3}' \
    > ${outputdir}/${sampleid}.reverse_Nucleosome.bed

  # Sort your generated BED files
  sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.forward_Nucleosome.bed -o ${outputdir}/${sampleid}.sortedNuc.forward_Nucleosome.bed
  sort -k 1V,1 -k 2n,2 ${outputdir}/${sampleid}.reverse_Nucleosome.bed -o ${outputdir}/${sampleid}.sortedNuc.reverse_Nucleosome.bed

  ##### sort with -k 1V,1 to get the correct order of chromosome, add -t first to get first nucleosome only, match row numbers. 
  # option "-t first" or "-t last". default "-t all". First intersected or last intersected region, or all regions. 
  bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.forward_Nucleosome.bed -b ${nucleosome_ref} -t first | awk -v OFS='\t' '{$13=$12 - $2;print $0}' > ${outputdir}/${sampleid}.forward_Nucleosome.dist.bed
  bedtools closest -a ${outputdir}/${sampleid}.sortedNuc.reverse_Nucleosome.bed -b ${nucleosome_ref} -t first | awk -v OFS='\t' '{$13=$12 - $2;print $0}' > ${outputdir}/${sampleid}.reverse_Nucleosome.dist.bed

  echo -e "sorting forward nucleosome file"
  sort -k4,4 ${outputdir}/${sampleid}.forward_Nucleosome.dist.bed > ${outputdir}/${sampleid}.forward_Nucleosome.dist.sorted.bed
  echo -e "sorting reverse nucleosome file"
  sort -k4,4 ${outputdir}/${sampleid}.reverse_Nucleosome.dist.bed > ${outputdir}/${sampleid}.reverse_Nucleosome.dist.sorted.bed

  touch ${outputdir}/${sampleid}.finished_Nucleosome.txt
fi

count_nuc_forward=$(cat ${outputdir}/${sampleid}.forward_Nucleosome.dist.sorted.bed | wc -l)
count_nuc_reverse=$(cat ${outputdir}/${sampleid}.reverse_Nucleosome.dist.sorted.bed | wc -l)

#####----------------------------------------------------------------------#####
##### Merge all features into one single tsv output file. 
#####----------------------------------------------------------------------#####
if [ ! -f "${outputdir}/${sampleid}.frag_output.tsv" ]; then
    echo -e "Merge all intermediate files into one single tsv output file"

    ##### column $13: distance of forward read to the nearest nucleosome
    cat ${outputdir}/${sampleid}.forward_Nucleosome.dist.sorted.bed | cut -f13 > ${outputdir}/forward_nuc.tmp.txt

    ##### column $13: distance of reverse read to the nearest nucleosome
    cat ${outputdir}/${sampleid}.reverse_Nucleosome.dist.sorted.bed | cut -f13 > ${outputdir}/reverse_nuc.tmp.txt
    
    ##### column $2: 4bp end motif from forward reads
    cat ${outputdir}/${sampleid}.forward_endmotif4bp.sorted.txt | cut -f2 > ${outputdir}/forward_4bpEM.tmp.txt
    
    # ##### column $13: 4bp end motif from reverse reads
    cat ${outputdir}/${sampleid}.reverse_endmotif4bp.sorted.txt | cut -f2 > ${outputdir}/reverse_4bpEM.tmp.txt

    paste ${inputfrag} ${outputdir}/forward_nuc.tmp.txt ${outputdir}/reverse_nuc.tmp.txt ${outputdir}/forward_4bpEM.tmp.txt ${outputdir}/reverse_4bpEM.tmp.txt > ${outputdir}/${sampleid}.frag_output.tsv
else
    echo -e "The input fragment file " ${inputfrag} " has been converted to tsv output file. Skip this step ..."
fi

if [ "${cleanup}" = "true" ]; then
  rm -rf ${outputdir}/*.txt;
  files_to_delete=$(ls ${outputdir} | grep -v '.*\.tsv$' | grep -v 'full_Nucleosome.dist.final.bed');
  for file in $files_to_delete; do
    rm -rf ${outputdir}/${file};
  done
elif [ "${cleanup}" = "false" ]; then
  echo -e "keep all intermediate files"
else
    echo "Invalid value for boolean flag: ${cleanup}. Use 'true' or 'false' only."
    exit 1
fi

echo -e "EM-forward count" "\t" $count_4bpEM_forward >> ${outputdir}/${sampleid}.sanity_check.tsv;
echo -e "EM-reverse count" "\t" $count_4bpEM_reverse >> ${outputdir}/${sampleid}.sanity_check.tsv;
echo -e "Nucleosome-forward count" "\t" $count_nuc_forward >> ${outputdir}/${sampleid}.sanity_check.tsv;
echo -e "Nucleosome-reverse count" "\t" $count_nuc_reverse >> ${outputdir}/${sampleid}.sanity_check.tsv;