#! /usr/bin/env python
# Time-stamp: <08-06-2026 m.utrosa@bcbl.eu>
"""
Extract ROI time series from beta images and plot a single violin 
plot per ROI, where each beta is an average from all runs for that
subject.

Before running this script ensure that you have resampled the atlas
correctly to the resolution of the functional images. See $atlas_path.
"""

# 00. Prerequisites -------------------------------------------------
# Import python packages
import bids
import warnings
import numpy as np
import nibabel as nib
from pathlib import Path

# Import custom-made functions
import grabber
import roisExtVis as rem

# 01. Parameters ----------------------------------------------------
# Experiment info
subID  = 5
sesIDs = [2, 3, 4, 5, 6, 7] # 2, 3, 4, 5, 6, 7
anatID = 2
acqIDs = ["BLOCK4"] # "BLOCK2", "BLOCK3", "BLOCK4"
space = "T1w"
conditions = [4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125] 
save = False
plot_rois = ["IC-L", "IC-R", "MGB-L", "MGB-R"] 
denoising = True

# Project directories
homePath   = Path("/home/mutrosa/mutrosa/Documents/devLoc")
dataPath   = homePath / "results" / "timDev_abs" / f"NORDIC-{denoising}" / "1stLevel"
atlas_name = f"sub-invivo_resampled_to-{space}_desc-betaNORDIC{denoising}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
atlas_path = homePath / "templates" / atlas_name
out_dir    = homePath / "results" / "visualization"
out_dir.mkdir(exist_ok=True)

# ROIs --------------------------------------------------------------
# Size represents the volume (mm3) of subcortical structures based on
# in-vivo functional clusters (see Table 1, Sitek et al., 2019).
# Label is identified from plotting unique atlas values in freeview.
rois = {'CN-L'  : {'size': 11, 'label': 1},
		'CN-R'  : {'size': 11, 'label': 2},
		'SOC-L' : {'size': 29, 'label': 3},
		'SOC-R' : {'size': 29, 'label': 4},
		'IC-L'  : {'size': 146, 'label': 5},
		'IC-R'  : {'size': 146, 'label': 6},
		'MGB-L' : {'size': 152, 'label': 7},
		'MGB-R' : {'size': 152, 'label': 8}}

# 02. Extract data from ROIs ----------------------------------------
summed_contrasts = {name: {c: [] for c in conditions }  for name in rois.keys()}
beta_affine = None

for sesID in sesIDs:
	for acqID in acqIDs:

		# Inspect the SPM.xX.name to see which beta images correspond
		# to which conditions of the SPM design matrix. 
		# Assuming numerical naming of beta images:		
		# beta_space-T1w_0004.nii				
		for b in range(1, len(conditions) + 1): 
			
			# Current condition
			cond = conditions[b - 1]

			# Construct the path
			name_beta = f"beta_space-{space}_{b:04d}.nii"
			beta_fold = dataPath / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
			beta_path = beta_fold / name_beta
			print(beta_path)

			# Extract the data		
			mask, _, beta_affine = rem.extract_roi_array(
				subID,
				sesID,
				acqID,
				atlas_path,
				space,
				beta_path,
				rois,
				out_dir,
				verbose=False,
				save=False
			)
		
			# Accumulate arrays for summation
			for name in rois.keys():
				summed_contrasts[name][cond].append(mask[name])

# 03. Average -------------------------------------------------------
final_betas = {name: {c: [] for c in conditions }  for name in rois.keys()}

# Iterate through the ROIs
for roi_name in rois.keys():

	roi_dict = summed_contrasts[roi_name]

	# Iterate through timing deviancy conditions	
	for condition in conditions:

		# Get a list of arrays for current ROI	
		array_lists = roi_dict[condition]

		# Avoid "division by zero" error by removing empty arrays	
		valid_arrays = [arr for arr in array_lists if arr.size > 0]
		if len(valid_arrays) != len(array_lists):
			message = f"For ROI {roi_name}, removing {len(valid_arrays) - len(array_lists)} empty array(s) before averaging."
			warnings.warn(message)
		
		if len(valid_arrays) > 0:

			# Assuming that every array has the same length	
			# Average arrays per timing deviandy condition
			# Note: len(arrays) == len(sesIDs) * len(acqIDs) unless empty arrays are removed
			averaged_array = np.average(array_lists, axis=0)
			final_betas[roi_name][condition].append(averaged_array)

			# Optionally save all averaged betas to disk
			if save:
				summed_filename = f"sub-{subID:02d}_roi-{roi_name}_space-{space}_cond-{condition}_type-average.nii.gz"
				summed_path = Path(out_dir) / summed_filename
				
				# Save the summed array
				if beta_affine is not None:
					nib.save(nib.Nifti1Image(averaged_array, beta_affine), summed_path)

# 04. Plot ----------------------------------------------------------
rem.plot_violins_average(final_betas, subID, plot_rois, out_dir, space, scale=False, save=True)