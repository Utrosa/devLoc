#!/usr/bin/env bash
# Time-stamp: <12-06-2026 m.utrosa@bcbl.eu>

set -eo pipefail
# -e => exits if any of the processes called generate a non-zero return code at the end.
# -o pipefail => deals with failures in the middle of a pipeline.

# Run the code in an environment specific to the project
source activate localizer_fMRI

# Subject-specific parameters
subID=5
anatID=2
project="devLoc"
n_noise_scans=1
nordic=1 # If 1, applies NORDIC denoising on functional BOLD images. If 0, skips NORDIC.

# Task-specific parameters: task and acqIDs must correspond
#task="localizer"
#acqIDs=("FUNCLOC")

#task="devLoc" # for non-split physio (sessions 2 & 6)
task="timDev"
#task="freqDev"
#acqIDs=("BLOCK1" "BLOCK2" "BLOCK3" "BLOCK4")

# To which acqIDs do we apply TAPAS?
# Tapas will find all funcional scans in the dicom folder and
# will list them based on the acquisition label in alphabetical order.
acqIDXs=("1" "2" "3" "4") # "BLOCK1" "BLOCK2" "BLOCK3" "BLOCK4" "FUNCLOC"

homePath="/home/mutrosa/mutrosa/Documents/devLoc"
funcPath="$homePath/data_MRI/sourcedata/denoised"

# Loop through the sessions
sessions=("3" "4" "5" "7")
for sesID in "${sessions[@]}"; do

	# STEP 0
	# Generate sidecar files to set up the configuration file.
	# Run script pre_import.py in the terminal.

	# STEP 1
	# a.) BIDSifies sourcedata (dicoms).
	# b.) Removes background noise from MP2RAGE UNI images (T1w).
	# c.) Applies denoising to the functional scans using NORDIC method.
	# d.) Removes the noise scan from functional scans in /func folder.

	# import_MRI:
	# If --bidsify flag is not included, the dicoms are not bidsified.
	# If --tidy flag is included, original BOLD images are  overwritten with nordic-denoised outputs.

	#echo "**************** STEP 1: Starting curation of MRI data *****************"

	#if [[ "$nordic" -eq 1 ]]; then
	#	python -m scripts.import.import_MRI \
	#			"$subID" "$sesID" "$anatID" "$task" "$project" \
	#			"$homePath" "${acqIDs[@]}" \
	#			--bidsify \
	#			--nordic

	#elif [[ "$nordic" -eq 0 ]]; then
	#	python -m scripts.import.import_MRI \
	#			"$subID" "$sesID" "$anatID" "$task" "$project" \
	#			"$homePath" "${acqIDs[@]}" \
	#			--bidsify
	#fi

	#python -m scripts.import.remove_noise_scan "$homePath" "$funcPath" "$subID" "$sesID" "$task" "$n_noise_scans" --overwrite

	#echo "*************************** Completed STEP 1 ***************************"

	# STEP 2: EVENTS
	# a.) Renames logfiles according to BIDS with info about the acquisition ID.
	#	  CRITICAL: acquisition labels MUST BE in the order of data collection !
	# b.) Copies behavioral logfiles to the "raw" BIDS-compliant folder in BIDS format.

	#echo "*********************** STEP 2: Moving LOGFILES ************************"

	# python -m scripts.import.import_LOG \
	#		"$homePath" "$subID" "$sesID" "$task" "${acqIDs[@]}"

	#echo "*************************** Completed STEP 2 ***************************"

	# STEP 3: BIOPAC
	# a.) Converts sourcedata (.acq) into TAPAS-compatible data (.mat or .txt).
	# b.) Preprocesses the compatible data.
	# c.) Calculates regressors.

	echo "******************** STEP 3: Starting PHYSIO import ********************"

	python -m scripts.import.import_PHYSIO "$subID" "$sesID" "$project" "$task" "$homePath" "${acqIDXs[@]}"
done
echo "*************************** Completed STEP 3 ***************************"

conda deactivate
