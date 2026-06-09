#! /usr/bin/env python
# Time-stamp: <05-06-2026 m.utrosa@bcbl.eu>
# -----------------------------------------------------------------------------
# Extracts and plots values from collected data within ROI masks from Sitek's 
# atlas. DOI: 10.7554/eLife.48932
# -----------------------------------------------------------------------------

# Import python packages
import bids
import pandas as pd
import numpy as np
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Import custom-made functions
import grabber

def extract_roi_array(subID, sesID, acqID, atlas, space, res_path, rois, out_dir, verbose, save):
	'''
	Extracts result values from the specified regions of interest (ROIs).
	All extracted ROIs are saved in out_dir.

	Parameters:
	- subID: integer number, identifying the participant
	- sesID: integer number, identifying the session info
	- acqID: string, identifying the functional MRI sequence
	- atlas: string, path to an established atlas
	- space: coordinate space of the input and output data (native T1w or MNI)
	- res_path: string, path to outputs of 1st Level Analysis with SPM in Nipype
	- rois: dictionary, specifying names, volume and atlas label of target ROIs.
	- out_dir: string, specifying the folder name for saving the results as .nii.gz
	- verbose: If True, prints affines and shape of atlas and result data in the terminal.
	- save: If True, saves the extracted roi array as a nifti file to disk.
 
	Returns:
	- res_rois: dictionary, extracted result values per each ROI
	- res_roi_paths: dictionary, paths to the extracted result values per each ROI

	'''
	# Load atlas image: Sitek
	atlas_img    = nib.load(atlas)
	atlas_data   = atlas_img.get_fdata()
	atlas_affine = atlas_img.affine

	# Load image from the analysis: result
	res_img    = nib.load(res_path)
	res_data   = res_img.get_fdata() 
	res_affine = res_img.affine

	
	# Compare affines and ignore tiny differences due to floating points
	if verbose:
		print("atlas shape\n", atlas_img.shape)
		print("\ninput shape\n", res_img.shape)
		print("\natlas affine\n",   atlas_affine)
		print("\ninput affine\n", res_affine)
		print("\n\nqform res\n", res_img.header.get_qform()[0])
		print("\nqform Sitek\n",  atlas_img.header.get_qform()[0])
		print("\n\nsform res\n", res_img.header.get_sform()[0])
		print("\nsform Sitek\n",  atlas_img.header.get_sform()[0])

	# Extract values per ROI
	res_rois = {}
	res_roi_paths = {}
	for name, roi in rois.items():
		mask_array  = (atlas_data == roi['label']).astype(float)
		contrast_array = res_data[mask_array > 0].flatten() # Alejandro says: "Plot this!"

		# Save result as zipped nifti
		if save:
			res_masked = mask_array * res_data
			result_filename = f"sub-{subID:02d}_ses-{sesID:02d}_acq-{acqID}_roi-{name}_space-{space}.nii.gz"
			result_path = out_dir / result_filename
			res_roi_paths[name] = result_path
			nib.save(nib.Nifti1Image(res_masked, res_affine), result_path)

		# Save voxels values for plotting
		res_rois[name] = contrast_array

	return res_rois, res_roi_paths, res_affine

def plot_violins(mask_paths, subID, sesID, acqIDs, out_dir, space, scale):

	rows = []
	for acq_name, roi_masks in mask_paths.items():
		
		for roi_name, roi_path in roi_masks.items():
			mask_img = nib.load(roi_path)
			vals = mask_img.get_fdata().flatten()
			vals = vals[vals != 0]
			if len(vals) < 5:
				print(f"Warning: very few voxels for {roi_name}, {acq_name}")
			rows.extend([{"ROI": roi_name, "acqID": acq_name, "values": v} for v in vals])
		
	df = pd.DataFrame(rows)
	print(df.head())

	for roi, group in df.groupby("ROI"):
		n_acq = len(acqIDs)
		fig, axes = plt.subplots(1, n_acq, figsize = (1.5 * n_acq, 8), sharey = True)
		
		if n_acq == 1:
			axes = [axes]
		
		for ax, acq in zip(axes, acqIDs):
			sub_df = group[group["acqID"] == acq]
			color_map = dict(zip(acqIDs, sns.color_palette("pastel", n_colors=len(acqIDs))))
			if not sub_df.empty:
				sns.violinplot(
					y = "values",
					data = sub_df,
					ax = ax,
					hue="ROI",
					legend = False,
					inner = "point",
					cut = 0, 
					palette = [color_map[acq]],
					bw_adjust = 0.5
				)

				ax.set_title(f"{acq}", fontsize = 8)
				ax.set_xlabel("")
			ax.set_xticks([])
			if scale == True:
				ax.set_ylim(-5, 10)

		fig.suptitle(f"sub-{subID:02d}_ses-{sesID:02d}_roi-{roi}", fontsize = 12)
		fig.tight_layout()
		fig_name = f"sub-{subID:02d}_ses-{sesID:02d}_roi-{roi}_space-{space}_violins.png"
		fig_path = out_dir / fig_name
		plt.savefig(fig_path, dpi = 200, bbox_inches = "tight")
		plt.close(fig)

def plot_violins_average(betas, subID, plot_rois, out_dir, space, scale, save):
	"""
	Extracts result values from the specified regions of interest (ROIs).
	All extracted ROIs are saved in out_dir.

	Parameters:
	- betas: a nested dict with beta arrays per timing deviancy and ROIs
	- subID: integer number, identifying the participant
	- plot_rois: a list of ROI labels. These will be the plotted ROIs.
	- out_dir: string, specifying the folder name for saving the results as .nii.gz
	- space: coordinate space of the input and output data (native T1w or MNI)
	- scale: If True, all subplots share the same y axis (scaled).
	- save: If True, saves the figure to disk.

	Returns:
	- Figure with violin subplots.
	"""

	# Initialize a list to store values	
	rows = []

	# Iterate through each roi
	for roi in betas.keys():

		# Only plot data for the selected regions
		if roi in plot_rois:

			# Get data per timing deviancy condition
			for timDev, beta_array in betas[roi].items():
				vals = beta_array[0]
				rows.extend([{"ROI": roi, "timDev": timDev, "values": v} for v in vals.flatten()])

		# Create a dataframe suitable for plotting		
		df = pd.DataFrame(rows)
		print(df.head())

	# Get unique time values
	rois = df["ROI"].unique()
	timDevs = df["timDev"].unique()
	violins = sns.color_palette("Set2", n_colors=len(timDevs))
	subplots = sns.color_palette("Set2", n_colors=len(rois))

	# Create a grid of subplots
	n_cond = int(len(plot_rois) / 2) # TODO: adapt for uneven numbers
	fig, axes = plt.subplots(n_cond, n_cond, figsize=(12, 10), sharey=True)
	axes = axes.flatten()

	# Plotting
	for i, roi in enumerate(plot_rois):
		ax = axes[i]
		
		# Filter data for the current ROI
		roi_data = df[df["ROI"] == roi]
		sns.violinplot(
			x = "timDev",
			y = "values", 
			data = roi_data,
			ax = ax,
			hue = "timDev",
			legend = False,
			inner = "stick", # Show individual observations
			cut = 0, 
			palette = violins,
			bw_adjust = 0.5
		)

		ax.set_title(f"{roi}", fontsize = 12)
		ax.set_xlabel("Timing Deviation [msec]", fontsize = 12)
		ax.set_ylabel("Mean Signal Change [β]", fontsize = 12)

		if scale == True:
			ax.set_ylim(-10, 10)

		# Hide any unused subplots if there are fewer than 4 ROIs
		for j in range(len(plot_rois), len(axes)):
			fig.delaxes(axes[j])


		fig.suptitle(f"sub-{subID:02d}", fontsize = 16, fontweight = "bold" )
		fig.tight_layout()

	# Optionally save
	if save:
		fig_name = f"sub-{subID:02d}_space-{space}_violins.png"
		fig_path = out_dir / fig_name
		plt.savefig(fig_path, dpi = 300, bbox_inches = "tight")
		plt.close(fig)

	# Show the plot
	plt.show()

# TODO: plot which plots IC-R and MGB-R on the same plot 