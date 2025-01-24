#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --mem=60GB

#note - have to export "label" variable

pipelinedir=/home/gmh4/projects/binding_affinity_md_af/af3_protein_rna_pipeline
output_dir=output
json_path=${output_dir}/${label}/${label}.input.json

#This step writes to the following file due to af3 logic
output_json=${output_dir}/${label}/${label}_data.json

if [ ! -e ${output_json} ];then
    $pipelinedir/run-af3-structure-prediction.sh \
    --json_path=${json_path} \
    --output_dir=$output_dir \
    --run_data_pipeline=true \
    --run_inference=false 
fi 

job_name=af3_${label}
out_log=$(dirname $json_path)/${job_name}.af3.%j.log

if [ -e ${output_json} ]; then
  sbatch --job-name $job_name -o $out_log --export json_path=$output_json run-af3-inference.sbatch
fi
