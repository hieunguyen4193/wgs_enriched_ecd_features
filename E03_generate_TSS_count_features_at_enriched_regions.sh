##### clone github pandepth: https://github.com/HuiyangYu/PanDepth.git
# git clone https://github.com/HuiyangYu/PanDepth.git 
while getopts "i:b:o:p:" opt; do
    case $opt in
        i) inputbam="$OPTARG"
        ;;
        b) inputbed="$OPTARG"
        ;;
        o) outputdir="$OPTARG"
        ;;
        p) pandepth="$OPTARG"
        ;;
        r) run_name="$OPTARG"
        ;;
        \?) echo "Usage: $0 -b <inputbam> -d <inputbed> -o <outputdir> -p <pandepth>"
                exit 1
        ;;
    esac
done

if [[ -z "$inputbam" || -z "$inputbed" || -z "$outputdir" || -z "$pandepth" ]]; then
    echo "Missing required arguments."
    echo "Usage: $0 -b <inputbam> -d <inputbed> -o <outputdir> -p <pandepth>"
    exit 1
fi

echo -e "working on sample " $inputbam; 
echo -e "using bed file " $inputbed;
echo -e "output will be saved to " $outputdir;

mkdir -p ${outputdir}/${run_name};

filename=$(echo $inputbam | xargs -n 1 basename);
filename=${filename%.bam*};

${pandepth} -i ${inputbam} -b ${inputbed} -o ${outputdir}/${run_name}/${filename}.bed;

# EOF
