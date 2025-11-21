##### clone github pandepth: https://github.com/HuiyangYu/PanDepth.git
# git clone https://github.com/HuiyangYu/PanDepth.git 
while getopts "b:d:o:p:" opt; do
    case $opt in
        b) input_bam="$OPTARG"
        ;;
        d) input_bed="$OPTARG"
        ;;
        o) output_dir="$OPTARG"
        ;;
        p) pandepth="$OPTARG"
        ;;
        \?) echo "Usage: $0 -b <input_bam> -d <input_bed> -o <output_dir> -p <pandepth>"
                exit 1
        ;;
    esac
done

if [[ -z "$input_bam" || -z "$input_bed" || -z "$output_dir" || -z "$pandepth" ]]; then
    echo "Missing required arguments."
    echo "Usage: $0 -b <input_bam> -d <input_bed> -o <output_dir> -p <pandepth>"
    exit 1
fi

echo -e "working on sample " $input_bam; 
echo -e "using bed file " $input_bed;
echo -e "output will be saved to " $output_dir;

mkdir -p ${output_dir};

filename=$(echo $input_bam | xargs -n 1 basename);
filename=${filename%.bam*};

${pandepth} -i ${input_bam} -b ${input_bed} -o ${output_dir}/${filename}.bed;