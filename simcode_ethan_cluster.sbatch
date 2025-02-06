#!/bin/bash
#
#SBATCH --array=1-2%2
#SBATCH --mem=5000
#SBATCH -t 04:00:00
#SBATCH -J mysim_test%A
#SBATCH -p secondary
#SBATCH -o mysim_test%A_%a.o%j
#SBATCH -e mysim_test%A_%a.e%j
#SBATCH --mail-user=emullen2@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# NOTE: %A = $SLURM_ARRAY_JOB_ID
# NOTE: %a = $SLURM_ARRAY_TASK_ID
./mysim --jobid "$SLURM_ARRAY_JOB_ID" --taskid "$SLURM_ARRAY_TASK_ID"