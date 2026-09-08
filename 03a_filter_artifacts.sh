#! /usr/bin/env bash
# Time-stamp: <08-09-2026 m.utrosa@bcbl.eu>
set -eo pipefail
# -e => exits if any of the processes called generate a non-zero return code at the end.
# -o pipefail => deals with failures in the middle of a pipeline.

# Run the code in an environment specific to the project.
# On Citrix, use `source activate`; elsewehere `conda activate`.
source activate localizer_fMRI

# Job-specific parameters
subID=5
acqIDs=("FUNCLOC") # "FUNCLOC" "BLOCK1" "BLOCK2" "BLOCK3" "BLOCK4"
task="localizer"   # timDev or localizer or devLoc
denoising="True"  # NORDIC applied or not during preproc
biopac=1 # exclude (0) or include (1) BIOPAC physiological regressors

# Project-specific directories
homePath='/home/mutrosa/mutrosa/Documents/projects/devLoc' # Citrix
#homePath='/home/mutrosa/Documents/projects/devLoc'           # Local
mriPath="$homePath/data_MRI/derivatives/NORDIC-$denoising/derivatives" # NORDIC-True or False
physioPath="$homePath/data_physio/raw/"

# -------------------------- Confound keys preferences  ---------------------------
# confound_keys=None # Will default to: trans, rot, csf, wm
# confound_keys=( \
#  			"csf" \
#  			"csf_derivative1" \
#  			"csf_derivative1_power2" \
#  			"csf_power2" \
#  			"white_matter" \
#  			"white_matter_derivative1" \
#  			"white_matter_derivative1_power2" \
#  			"white_matter_power2" \
#  			"csf_wm"
#  		)
# The rigid body keys must be in order in which FSL expects them
# https://fsl.fmrib.ox.ac.uk/fsl/docs/registration/mcflirt.html
confound_keys=( \
 	'rot_x' \
 	'rot_y' \
 	'rot_z' \
 	'trans_x' \
 	'trans_y' \
 	'trans_z'
 	)

# ------------------------------ Filter Artifacts ------------------------------ 
## a.) Extracts specified confounds and motion artifacts from fMRIprep derivatives.
## b.) Optionally adds physiological regressors (TAPAS) to the confounds dataframe.

echo "**************** STEP 1: Filtering confounds & artifacts ***************"
for sesID in 2 6 ; do
	if [[ "$biopac" -eq 1 ]]; then
		python -m scripts.analysis.filter_artifacts \
				"$homePath" "$mriPath" "$physioPath" \
				"$subID" "$sesID" "$task" "$denoising" \
				--acqIDs "${acqIDs[@]}" \
				--confound_keys "${confound_keys[@]}" \
				--include_biopac

	elif [[ "$biopac" -eq 0 ]]; then
		python -m scripts.analysis.filter_artifacts \
			"$homePath" "$mriPath" "$physioPath" \
			"$subID" "$sesID" "$task" "$denoising" \
			--acqIDs "${acqIDs[@]}" \
			--confound_keys "${confound_keys[@]}"
	fi
done
echo "*************************** Completed STEP 1 ***************************"

conda deactivate
