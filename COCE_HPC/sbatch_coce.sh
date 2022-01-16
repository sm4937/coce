#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=eowm
#SBATCH --mail-type=END
#SBATCH --output=output/output_%A_%a.out
#SBATCH --error=error/error_%A_%a.err


module purge
module load matlab/2020b

export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
export MATLAB_LOG_DIR=$SLURM_JOBTMP

cd /scratch/sm4937/COCE/HBI
/share/apps/matlab/2020b/bin/matlab -nodisplay -r "execute_HBI_HPC" >> Log/Log_${SLURM_ARRAY_TASK_ID}.txt
#/share/apps/matlab/2020b/bin/matlab -nodisplay -r "execute_TAFKAP_HPC(${SLURM_ARRAY_TASK_ID})" >> Log/Log_${SLURM_ARRAY_TASK_ID}.txt

#rm -r $SCRATCH/MATLABJobStorage/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
