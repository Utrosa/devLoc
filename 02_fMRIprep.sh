#!/bin/bash
# Time-stamp: <09-05-2026>

# ------------------------ SGE options ------------------------------
#$ -q long.q
#$ -cwd
#$ -N devLoc_sub5
#$ -m ea
#$ -M m.utrosa@bcbl.eu
#$ -o preproc_logs/devLoc_output.txt
#$ -e preproc_logs/devLoc_error.txt
#$ -S /bin/bash

# ------------------------ JOB execution ----------------------------
# Set monitoring for errors
set -euo pipefail

echo "***** Job started *****"
date

# Load modules and list them for reproducibility
module load singularity/3.7.0
module list

# Subjects to process
subjects=("05")

# Run the fMRIprep preprocessing pipeline
for subID in "${subjects[@]}"; do
	bash data_MRI/code/preproc_singleSub_singularity.sh "$subID"
done

echo "**** Job ended ****"
date
