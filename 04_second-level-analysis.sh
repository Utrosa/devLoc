#! /usr/bin/env bash
# Time-stamp: <21-08-2026 m.utrosa@bcbl.eu>
# 2nd level GLM analysis on contrast and beta images
set -eo pipefail

# Run the code in an environment specific to the project.
source activate nipypee

#----------------------------------------------------------------------------------

# STEP 1: CONTRASTS
# a.) Obtains contrast images.
# b.) Inferential statistics on the extracted contrasts.
# c.) Visualization of the stats output.

echo "******** STEP 1: SPM 2nd level GLM analysis on contrast images ********"

# python -m scripts.analysis.contrast_get \
# 	"$homePath" \
# 	"$MNIpath" \
# 	--subIDs "$subIDs" \
# 	--sesIDs "$sesIDs" \
# 	--acqIDs "${acqIDs[@]}" \
# 	--smoothing "$smoothing"

echo "*************************** Completed STEP 1 ***************************"

# STEP 2: BETA COEFFICIENTS
# a.) Obtains beta images.
# b.) Inferential statistics on the extracted contrasts.
# c.) Visualization of the stats output.

echo "******** STEP 1: SPM 2nd level GLM analysis on beta images ********"

# ts=$(python -m scripts.analysis.roi_extraction \
# 	"$subID" \
# 	"$sesID" \
# 	"$ts" \
# 	"$homePath" \
# 	"$atlasPath" \
# 	"${acqIDs[@]}")

echo "*************************** Completed STEP 2 ***************************"
