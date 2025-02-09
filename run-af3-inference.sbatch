#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=60GB
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|h100
#SBATCH --job-name=run-inference

set -x 
pipelinedir=/home/gmh4/projects/binding_affinity_md_af/af3_protein_rna_pipeline

#NOTE - must set json_path from outside
label=$(basename $json_path .data.json)
model_cif=$(dirname $json_path)/${label}_model.cif
out_dir=$(dirname $json_path)

if [ ! -e "${model_cif}" ];then

    $pipelinedir/run-af3-structure-prediction.sh \
    --json_path=${json_path} \
    --output_dir=$out_dir \
    --run_data_pipeline=false \
    --run_inference=true

fi

job_name=af3_md_${label}
out_log=$(dirname $json_path)/${job_name}.md.%j.log

md_dir=$(out_dir)/md_implicit
cp ${model_cif} $md_dir
md_cif=${md_dir}/${label}_model.cif

if [ -e ${md_cif} ]; then
  sbatch --job-name $job_name -o $out_log --export model_cif=$md_cif run-md-af-model.sbatch
fi
