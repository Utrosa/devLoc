#! /usr/bin/env python
# Time-stamp: <07-09-2026 m.utrosa@bcbl.eu>
"""
Extracting values from collected data within ROI masks from Sitek's atlas
Plot a single violin plot per ROI, where the betas are an average from all
runs.

Before running this script ensure that you have resampled the atlas
correctly to the resolution of the functional images (ref: atlas_path)!
"""
# CHECK: mask_paths[acqID] = mask_path
import config as c
import numpy as np
import nibabel as nib
from utils import extract_roi_array, plot_violins

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Extract data from ROIs
# Initialize a dictionary to save extracted values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extracted_spmT = {acqID: {s: [] for s in c.sesIDs} for acqID in c.acqIDs}
roi_names = list(c.rois.keys())
all_individual_paths = {acqID: {s: [] for s in c.sesIDs} for acqID in c.acqIDs}

for acqID in c.acqIDs:
	for sesID in c.sesIDs:

		# Construct the path
		spmT_fold = c.dataPath / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
		spmT_path = spmT_fold / c.spmT_filename + "0001.nii"

		# Extract the subcortical arrays			
		masks_subcor, mask_path_subcor, spmT_affine_subcor = extract_roi_array(
			c.subID,
			sesID, 
			acqID, 
			c.atlas_subcor_path, 
			c.space, 
			spmT_path, 
			c.rois_subcortical, 
			c.out_1st,
			verbose=False,
			save=c.save_roi,
			average_voxels=False # Keeping this false for consistency
		)

		# Extract the cortical arrays
		masks_cor,  mask_path_cor, spmT_affine_cor = extract_roi_array(
			c.subID,
			sesID, 
			acqID, 
			c.atlas_cor_path, 
			c.space, 
			spmT_path, 
			c.rois_cortical, 
			c.out_1st,
			verbose=False,
			save=c.save_roi,
			average_voxels=False # Keeping this false for consistency
		)

		# Accumulate paths
		mask_path_all = mask_path_cor | mask_path_subcor
		all_individual_paths[acqID][sesID] = mask_path_all

		# Assign affine
		if spmT_affine_subcor.all() == spmT_affine_cor.all():
			spmT_affine = spmT_affine_subcor

		# Accumulate subcortical arrays for summation
		masks_all = masks_cor | masks_subcor
		extracted_spmT[acqID][sesID].append(masks_all)
		for roi, mask_path_ind in mask_path_all.items():
			if roi in c.rois_subcortical.keys():
				nib.save(nib.Nifti1Image(masks_subcor[roi], spmT_affine), mask_path_ind)

		# Accumulate cortical arrays for summation
		for roi, mask_path_ind in mask_path_all.items():
			if roi in c.rois_cortical.keys():
				nib.save(nib.Nifti1Image(masks_cor[roi], spmT_affine), mask_path_ind)

# Print shape of the raw extracted data
print("\nAssuming all beta images have the same affine.")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Transform the extracted values.
# Average across runs or voxels: (n_runs, n_voxels)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final_spmT_rois = {}
final_summed_paths = {}
for acqID in c.acqIDs:
	for sesID in c.sesIDs:
		for name in c.rois.keys():
			array_lists = extracted_spmT[acqID][sesID][0]
			array_list  = array_lists[name]

			if len(array_list) > 0:

				# Sum across sessions
				final_spmT_rois[name] = np.sum(array_list, axis=0)

				# Save the summed array
				if spmT_affine is not None:
					summed_filename = f"sub-{c.subID:02d}_ses-{c.sessions}_acq-{c.blocks}_roi-{name}_space-{c.space}.nii.gz"
					summed_path = c.out_1st / summed_filename
					final_summed_paths[name] = summed_path
					nib.save(nib.Nifti1Image(final_spmT_rois[name], spmT_affine), summed_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03. Descriptive plotting per acquisition (grouping across sessions)
# TODO: summation across sessions but keeping acquisitions
# TODO: averaging voxels
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting non-summed, voxelwise data
plot_violins(
	all_individual_paths, 
	c.subID, 
	c.sessions, 
	c.acqIDs, 
	c.out_1st, 
	c.space, 
	scale=True
)

# TODO: IMPROVE: Plotting summed (collapsing sessions & acquisitons), voxelwise data
# summed_paths = {c.blocks: {c.sessions : final_summed_paths}}
# plot_violins(
# 	all_individual_paths, 
# 	c.subID, 
# 	c.sessions, 
# 	c.acqIDs, 
# 	c.out_1st, 
# 	c.space, 
# 	scale=True
# )