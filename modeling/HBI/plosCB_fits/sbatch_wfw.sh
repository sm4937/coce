#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=99:00:00
#SBATCH --mem=40GB
#SBATCH --job-name=wfw
#SBATCH --mail-type=END
#SBATCH --mail-user=sm4937@nyu.edu
#SBATCH --output=output/output_%A_%a.out
#SBATCH --error=error/error_%A_%a.err


module purge
module load matlab/2020b

export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
export MATLAB_LOG_DIR=$SLURM_JOBTMP

cd /scratch/sm4937/delta_tests/modeling/HBI/quick_delta_tests
/share/apps/matlab/2020b/bin/matlab -nodisplay -r "run_model_fitting_deltas" >> Log/Log.txt

#rm -r $SCRATCH/MATLABJobStorage/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
# send a {SLURM_ARRAY_TASK_ID} by running sbatch sbatch_eowm.sh -f 4 ?