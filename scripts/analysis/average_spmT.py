#! /usr/bin/env python
# Time-stamp: <02-06-2026 m.utrosa@bcbl.eu>
"""
Extracting values from collected data within ROI masks from Sitek's atlas
Plot a single violin plot per ROI, where the betas are an average from all
runs.

Before running this script ensure that you have resampled the atlas
correctly to the resolution of the functional images (ref: atlas_path)!
"""

# 00. Prerequisites -------------------------------------------------
# Import python packages
import bids
import numpy as np
import nibabel as nib
from pathlib import Path

# Import custom-made functions
import grabber
import roisExtVis as rem

# 01. Parameters ----------------------------------------------------
# Experiment info
subID  = 5
sesIDs = [2, 3, 4, 5, 6, 7]
anatID = 2
acqIDs = ["BLOCK1", "BLOCK2", "BLOCK3", "BLOCK4"]
space = "T1w"

# Project directories
homePath   = Path("/home/mutrosa/mutrosa/Documents/projects/devLoc")
atlas_path = homePath / "templates" / f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
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
summed_contrasts = {name: [] for name in rois.keys()}
all_individual_paths = {}	

for sesID in sesIDs:
	for acqID in acqIDs:

		# Construct the path
		name_spmT   = f"spmT_space-{space}_0001_trans_out.nii.gz"
		folder_spmT = "_acqID_{acqID}_anatID_{anatID}_sesID_{sesID}_subID_{subID}"
		spmT_path   = homePath / "results" / "1stLevel" / folder_spmT / name_spmT

		# Extract the data for current session			
		masks, mask_path = extract_roi_array(subID, sesID, acqID, atlas_path, space, spmT_path, rois, out_dir)
		# CHECK: mask_paths[acqID] = mask_path
		temp_img = nib.load(spmT_path)
		last_affine = temp_img.affine

		# Accumulate arrays for summation
		for name in rois.keys():
			summed_contrasts[name].append(masks[name])
			all_individual_paths[acqID] = mask_path

# 03. Average -------------------------------------------------------
final_spmT_rois = {}
final_summed_paths = {}
for name in rois.keys():
	array_list = summed_contrasts[name]
	if len(array_list) > 0:
		final_spmT_rois[name] = np.sum(array_list, axis=0)

	     # Create filename for the summed result
		summed_filename = f"sub-{subID:02d}_ses-sum_acq-sum_roi-{name}_space-{space}.nii.gz"
		summed_path = Path(out_dir) / summed_filename
		
		# Save the summed array
		if last_affine is not None:
			nib.save(nib.Nifti1Image(final_spmT_rois[name], last_affine), summed_path)
			final_summed_paths[name] = summed_path

# 04. Plot per session ----------------------------------------------
mask_paths = {}
mask_paths[acqIDs[0]] = final_summed_paths
sesID = 234567 # all the sessions
plot_violins(mask_paths, subID, sesID, acqIDs, out_dir, space, scale=True)

mask_paths = all_individual_paths
for sesID in sesIDs:
	plot_violins(mask_paths, subID, sesID, acqIDs, out_dir, space, scale=True)