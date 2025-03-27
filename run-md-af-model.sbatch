#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=20GB
#SBATCH --gres=gpu:1
#SBATCH --job-name=run-md-from-cif

pipelinedir=/home/gmh4/projects/binding_affinity_md_af/af3_protein_rna_pipeline

#set model_cif from outside

openmm=$pipelinedir/run-openmm.sh 
$openmm python $pipelinedir/run_openmm_setup_md.py $model_cif
