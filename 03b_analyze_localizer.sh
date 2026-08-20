#! /usr/bin/env bash
# Time-stamp: <12-05-2026 m.utrosa@bcbl.eu>

set -eo pipefail
# -e => exits if any of the processes called generate a non-zero return code at the end.
# -o pipefail => deals with failures in the middle of a pipeline.

# Run the code in an environment specific to the project.
source activate localizer_fMRI

MNIpath="$homePath/templates/tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz"
atlasPath="$homePath/templates/atlas/invivo_resampled_to-MNI_res-01.nii.gz"
smoothing=2  # If larger than zero, applies smoothing kernel of that size in mm.
volterra=0   # 0: no volterra; 1: applies volterra

#----------------------------------------------------------------------------------

# STEP 1: Nipype 1st level GLM analysis
# a.) Specifies nodes, connects them, and visualizes the workflow.
# b.) Generates the SPM-specific model, the design matrix, and estimates contrasts.
# c.) Warps results from T1w space to MNI.

echo "***************** STEP 1: SPM 1st level GLM analysis *******************"

# python -m scripts.analysis.first_level_analysis \
# 	"$homePath" \
# 	"$MNIpath" \
# 	--subIDs "$subIDs" \
# 	--sesIDs "$sesIDs" \
# 	--acqIDs "${acqIDs[@]}" \
# 	--smoothing "$smoothing"

echo "*************************** Completed STEP 1 ***************************"

# STEP 2: Extraction of target ROIs & visualization
# a.) Extracts target ROIs
# b.) Visualizes

echo "***************** STEP 2: Extracting ROIs and plotting *****************"

# ts=$(python -m scripts.analysis.roi_extraction \
# 	"$subID" \
# 	"$sesID" \
# 	"$ts" \
# 	"$homePath" \
# 	"$atlasPath" \
# 	"${acqIDs[@]}")

echo "*************************** Completed STEP 2 ***************************"

