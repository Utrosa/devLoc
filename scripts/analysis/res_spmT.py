#! /usr/bin/env python
# Time-stamp: <31-08-2026 m.utrosa@bcbl.eu>
"""
Extracting values from collected data within ROI masks from Sitek's atlas
Plot a single violin plot per ROI, where the betas are an average from all
runs.

Before running this script ensure that you have resampled the atlas
correctly to the resolution of the functional images (ref: atlas_path)!
"""
# WHY IS THIS NOT DONE FOR CORTICAL ROIS?
# 00. Configuration -----------------------------------------------------------
# Import custom-made scripts
import grabber
import config as c
from utils import extract_roi_array, plot_violins

# 01. Extract data from ROIs --------------------------------------------------
summed_contrasts = {name: [] for name in c.rois_subcortical.keys()}
all_individual_paths = {}	

for sesID in c.sesIDs:
	for acqID in c.acqIDs:

		# Construct the path
		name_spmT   = f"spmT_space-{c.space}_0001_trans_out.nii.gz"
		folder_spmT = "_acqID_{acqID}_anatID_{c.anatID}_sesID_{sesID}_subID_{c.subID}"
		spmT_path   = homePath / "results" / "1stLevel" / folder_spmT / name_spmT

		# Extract the data for current session			
		masks, mask_path = extract_roi_array(c.subID, sesID, acqID, atlas_path, c.space, spmT_path, c.rois_subcortical, c.out_dir_spmt)
		# CHECK: mask_paths[acqID] = mask_path
		temp_img = nib.load(spmT_path)
		last_affine = temp_img.affine

		# Accumulate arrays for summation
		for name in c.rois_subcortical.keys():
			summed_contrasts[name].append(masks[name])
			all_individual_paths[acqID] = mask_path

# 02. Average -----------------------------------------------------------------
final_spmT_rois = {}
final_summed_paths = {}
for name in rois_subcortical.keys():
	array_list = summed_contrasts[name]
	if len(array_list) > 0:
		final_spmT_rois[name] = np.sum(array_list, axis=0)

	     # Create filename for the summed result
		summed_filename = f"sub-{subID:02d}_ses-sum_acq-sum_roi-{name}_space-{c.space}.nii.gz"
		summed_path = Path(c.out_dir_spmt) / summed_filename
		
		# Save the summed array
		if last_affine is not None:
			nib.save(nib.Nifti1Image(final_spmT_rois[name], last_affine), summed_path)
			final_summed_paths[name] = summed_path

# 04. Plot per session --------------------------------------------------------
mask_paths = {}
mask_paths[c.acqIDs[0]] = final_summed_paths
sesID = 234567 # all the sessions # TODO: this seems not too good ...
plot_violins(
	mask_paths, 
	c.subID, 
	sesID, 
	c.acqIDs, 
	c.out_dir_spmt, 
	c.space, 
	scale=True
)

mask_paths = all_individual_paths
for sesID in c.sesIDs:
	plot_violins(
		mask_paths, 
		c.subID, 
		sesID, 
		c.acqIDs, 
		c.out_dir_spmt, 
		c.space, 
		scale=True
	)