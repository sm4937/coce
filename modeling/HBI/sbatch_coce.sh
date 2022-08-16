#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --mem=100GB
#SBATCH --job-name=eowm
#SBATCH --mail-type=END
#SBATCH --mail-user=sm4937@nyu.edu
#SBATCH --output=output/output_%A_%a.out
#SBATCH --error=error/error_%A_%a.err


module purge
module load matlab/2020b

export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
export MATLAB_LOG_DIR=$SLURM_JOBTMP

cd /scratch/sm4937/coce/HBI
#/share/apps/matlab/2020b/bin/matlab -nodisplay -r "execute_TAFKAP_subj${SLURM_ARRAY_TASK_ID}" >> Log/Log_${SLURM_ARRAY_TASK_ID}.txt
# I think the way above is how to send an input from the command line
# to hack it for now I'm doing it this way
/share/apps/matlab/2020b/bin/matlab -nodisplay -r "HPC_runModelSearch_withNull" >> Log/Log_${ID}.txt

#rm -r $SCRATCH/MATLABJobStorage/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
# send a {SLURM_ARRAY_TASK_ID} by running sbatch sbatch_eowm.sh -f 4 ?