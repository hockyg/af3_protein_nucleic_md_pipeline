#!/bin/bash

args=
for i in "$@"; do 
  i="${i//\\/\\\\}"
  args="${args} \"${i//\"/\\\"}\""
done

cmd="python /ext3/alphafold3/run_alphafold.py"

if [[ "${args}" == "" ]]; then 
  cmd=
  args="/bin/bash"
fi

if [[ -e /dev/nvidia0 ]]; then nv="--nv"; fi

if [[ "${SINGULARITY_CONTAINER}" != "" ]]; then 
  export PATH=/share/apps/apptainer/bin:${PATH}
fi

ALPHAFOLDE3_DATA_HOME_host="/vast/work/public/alphafold-3.0.0-20241121"
mmcif_files_sqf="/vast/work/public/alphafold-3.0.0-20241121/public_databases/pdb_2022_09_28_mmcif_files.sqf"

ALPHAFOLDE3_DATA_HOME="/alphafolder3_data"

if [[ "${AF3_MODELS_DIR}" != "" ]] && [[ -d ${AF3_MODELS_DIR} ]]; then
  af3_models_bind="--bind ${AF3_MODELS_DIR}:${ALPHAFOLDE3_DATA_HOME}/models:ro"
fi

singularity exec ${nv} \
--overlay /scratch/work/public/apps/alphafold/3.0.0/alphafold-3.0.0-20241122.sqf:ro \
--bind ${ALPHAFOLDE3_DATA_HOME_host}:${ALPHAFOLDE3_DATA_HOME}:ro \
${af3_models_bind} \
--mount type=bind,src=${mmcif_files_sqf},dst=${ALPHAFOLDE3_DATA_HOME}/public_databases/mmcif_files,ro,image-src=/mmcif_files \
/scratch/work/public/singularity/cuda12.6.2-cudnn9.5.0-devel-ubuntu24.04.1.sif \
/bin/bash -c "
unset -f which
if [[ -e /ext3/env.sh ]]; then source /ext3/env.sh; fi
${cmd} ${args}
"
