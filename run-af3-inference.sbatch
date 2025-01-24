#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=60GB
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100|h100
#SBATCH --job-name=run-inference

pipelinedir=/home/gmh4/projects/binding_affinity_md_af/af3_protein_rna_pipeline

#NOTE - must set json_path from outside

$pipelinedir/run-af3-structure-prediction.sh \
--json_path=${json_path} \
--output_dir=output \
--run_data_pipeline=false \
--run_inference=true
