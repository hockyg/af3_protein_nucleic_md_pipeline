#!/bin/bash
input_file=$1
if [ -z "$input_file" ];then
	echo "Usage: $0 input_json_file"
        exit 1
fi

label=$(basename $input_file .input.json)

#Step 2, run alphafold prediction in two steps
job_name=predict_$label
sbatch --export input_json=$input_file --job-name $job_name -o ${prefix}.pipeline.%j.out run-af3-pipeline-full-fromjson.sbatch
