#!/bin/bash
input_file=$1
if [ -z "$input_file" ];then
	echo "Usage: $0 input_file"
	echo -e "\tFirst line = Protein sequence\n\tSecond line = DNA/RNA sequence"
        exit 1
fi

label=$(basename $input_file .input)
prefix=output/$label/$label
output_json=${prefix}.input.json

#Step 1, create output json from input file
if [ ! -e "$output_json" ];then
    python input_to_json.py $input_file $output_json
fi

#Step 2, run alphafold prediction in two steps
job_name=predict_$label
sbatch --export label=$label --job-name $job_name -o ${prefix}.pipeline.%j.out run-af3-pipeline-full.sbatch
