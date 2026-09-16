#! /usr/bin/env python
# Time-stamp: <15-09-2026 m.utrosa@bcbl.eu>
"""
Extract ROI arrays from beta images and plot a single violin 
plot per ROI, where each beta is:
	a.) an average across runs for that subject (n = n_voxels), or
	b.) an average across voxels per run (n = n_runs).

Before running this script ensure that you have resampled the atlas
and outputs from 1st level GLM analysis correctly to the desired space.
"""
# Import python packages
import pandas as pd
import numpy as np
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from itertools import combinations

# Import custom-made functions
import config as c
from config import plotConf, apply_figure_style
from utils import extract_roi_array, plot_violins_average

# Apply style for the figures
apply_figure_style()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 00. Check that inputs are defined correctly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Averaging configuration
if c.average_voxels and c.average_runs:
    raise ValueError(
        "Statistical tests cannot be performed when averaging across BOTH runs and voxels. "
        "Please set either average_voxels=False or average_runs=False."
    )

if not c.average_voxels and not c.average_runs:
    raise ValueError(
        "To perform statistical tests we need one-dimensional arrays, "
        "which means that the selected values (per beta image) "
        "have to be averaged EITHER across runs or voxels."
    )

# ROI specification (names and order)
roi_names = list(c.rois.keys())
print(f"The configured ROIS are: {roi_names}.")

# Selection of conditions: the order matters!
conditions = c.timDevs
print(f"\nThe selected {len(conditions)} conditions are:\n{conditions}."
	"\nIMPORTANT: The above conditions must correspond to the order of regressors"
	" in the SPM design. The order impacts which beta files are loaded.")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Extract beta values per voxel from each ROI.
# {ROI1: {C1: [[v1, v2, v3, ...],[v1, v2, v3, ...]], 
#         C2: [[v1, v2, v3, ...],[v1, v2, v3, ...]]} 
# Each roi is a dictionary of length n_cond
# Each cond is a list with n_run arrays
# Each run array has n_voxel values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize a dictionary to save extracted values
extracted_beta = {name: {c: [] for c in conditions} for name in roi_names}
beta_affine = None
for b in range(1, len(conditions) + 1):

    # Current condition
    cond = conditions[b - 1]

    # Construct the name
    beta_name = f"beta_{c.resampled_stem}_{b:04d}.nii"

    # Define beta paths
    # TODO: this doesn't care about sessions ...clear
    beta_paths = list(c.dataPath.rglob(beta_name))

    for beta_path in beta_paths:

        # Extract the subcortical arrays        
        mask_subcor, _, beta_subcor_affine = extract_roi_array(
            c.atlas_subcor_path,
            c.space,
            beta_path,
            c.rois_subcortical,
            c.out_1st,
            verbose=c.verbose,
            save=c.save_roi,
            average_voxels=False # Keeping this false for consistency
        )
        
        # Extract the cortical arrays
        mask_cor, _, beta_cor_affine = extract_roi_array(
            c.atlas_cor_path,
            c.space,
            beta_path,
            c.rois_cortical,
            c.out_1st,
            verbose=c.verbose,
            save=c.save_roi,
            average_voxels=False # Keeping this false for consistency
        )

        # Accumulate subcortical arrays for summation
        for name in c.rois_subcortical.keys():
            extracted_beta[name][cond].append(mask_subcor[name])

        # Accumulate cortical arrays for summation
        for name in c.rois_cortical.keys():
            extracted_beta[name][cond].append(mask_cor[name])

# Assign beta affine
if beta_cor_affine.all() == beta_subcor_affine.all():
	beta_affine = beta_subcor_affine
	print("\nAssuming all beta images have the same affine.")

# Print shape of the raw extracted data
print("\n--- RAW EXTRACTED DATA ---")
for roi_name in roi_names:
	for cond in conditions:
		print(f"{roi_name} - {cond}: {np.shape(extracted_beta[roi_name][cond])}")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Transform the extracted beta values.
# Average across runs or voxels. Original shape: (n_runs, n_voxels).
# If we have a single dataset (one subject, one session) with concatenation per
# block it will be (n_observations, n_voxels) per ROI.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
selected_betas = {name: {c: [] for c in conditions } for name in roi_names}

# Iterate through the ROIs
for roi_name in roi_names:
	roi_dict = extracted_beta[roi_name]

	# Iterate through the conditions (the SPM regressors)
	for condition in conditions:

		# Get a list of arrays (n_observations length)
		array_list = roi_dict[condition]

		# Optionally, remove empty arrays (zero values)
		if c.remove_empty:
			valid_array = [arr for arr in array_list if arr.size > 0]
			if len(valid_array) != len(array_list):
				message = f"\nFor ROI {roi_name}, removing {len(array_list) - len(valid_array)} empty array(s) before averaging."
				warnings.warn(message)
		else:
			valid_array = array_list

		if len(valid_array) > 0:

			# Collapse voxels: get a mean contrast value across voxels
			if c.average_voxels:
				averaged_array = np.mean(valid_array, axis=-1)
				selected_betas[roi_name][condition].append(averaged_array)

			# Collapse runs: get a mean contrast value across runs
			elif c.average_runs:
				if np.ndim(valid_array) == 1:
					selected_betas[roi_name][condition].append(valid_array)
				elif np.ndim(valid_array) == 2:
					averaged_array = np.mean(valid_array, axis=0)
					selected_betas[roi_name][condition].append(averaged_array)
				else:
					raise ValueError(
						f"Unusual number of dimensions for beta arrays: {np.ndim(valid_array)}"
						"Check the output data: how did you combine it over sessions, runs, subjects?")
			
			# Optionally save all averaged betas to disk
			# TODO: Why here averaged and the contrasts summed?
			if c.save_averaged:
				summed_filename = f"betas_roi-{roi_name}_sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_job-{c.jobName}_space-{c.resampled_stem}_cond-{condition}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}.nii.gz"
				summed_path = c.out_1st / summed_filename
				
				# Save the summed array
				if beta_affine is not None:
					nib.save(nib.Nifti1Image(averaged_array, beta_affine), summed_path)

# Print update on the structure of array
print("\n--- TRANSFORMED EXTRACTED DATA ---",
     f"\nAveraged across runs: {c.average_runs}",
     f"\nAveraged across voxels: {c.average_voxels}\n")
for roi_name in roi_names:
	for cond in conditions:
	    if np.isscalar(selected_betas[roi_name]):
	        print(f"{roi_name}, {cond}: {float(selected_betas[roi_name][cond]):.4f}")
	    else:
	        print(f"{roi_name}, {cond}: {np.shape(np.array(selected_betas[roi_name][cond]))}")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03. Descriptive plotting
# Plot beta values per ROI (subplot) and condition (x axis)
# AIM: visualize the range and shape of betas' distributions.
# TODO: assess significance, no? DIFFERENT FROM ZERO OR NOT? 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_violins_average(
		selected_betas,
		c.subID,
		c.plot_rois,
		c.plotConf,
		c.out_1st,
		c.space,
		scale=False,
		save=c.save_fig,
		show=c.show_fig,
		average_runs=c.average_runs,    # only for the filename
		average_voxels=c.average_voxels # only for the filename
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04a. Inferential statistics
# Compare conditions within an ROI (either across runs or voxels).
# RQ: Is there a significant difference between conditions in the ROI?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reshape the data for plotting
data_rows = []
for roi_name, cond_dict in selected_betas.items():
	for cond in conditions:
		values = cond_dict[cond][0]
		
		# Check dimensions in case of a single voxel/run
		if np.ndim(values) == 0:
			values_to_iterate = [values]
		else:
			values_to_iterate = values

		# Append info for the dataframe
		for idx, val in enumerate(values_to_iterate):
			data_rows.append({
				'ROI': roi_name,
				'Condition': cond,
				'Value': val,
				'Idx': idx + 1,
				'nVox_nRun': extracted_beta[roi_name][cond][0].shape,
				'avgVox' : c.average_voxels,
				'avgRun': c.average_runs
			})

# Create dataframe
df = pd.DataFrame(data_rows)

# Rename index column to reflect the condition that was not collapsed.
if c.average_runs:
	df_renamed = df.rename(columns={"Idx": "Voxel"}, inplace=False)
elif c.average_voxels:
	df_renamed = df.rename(columns={"Idx": "Run"}, inplace=False)
print(f"The renamed dataframe:\n{df_renamed.head()}")

# Perform the statistical tests for all pairs of conditions in the ROI
results_table = []
condition_pairs = list(combinations(conditions, 2))

# Pivot the data by IDX (run or voxel) and ROI
pivot_df = df.pivot_table(index=['ROI', 'Idx'], columns='Condition', values='Value')

# Iterate through the ROIs
for roi in roi_names:
	roi_data = pivot_df.loc[roi]

	for cond_a, cond_b in condition_pairs:
		pair_data = roi_data[[cond_a, cond_b]]

		# Optionally, remove nan values
		if c.remove_empty:
			valid_pair_data = pair_data.dropna()
		else:
			valid_pair_data = pair_data

		# Number of observations
		n_obs = len(valid_pair_data)

		# Extract the arrays
		values_a = valid_pair_data[cond_a].values
		values_b = valid_pair_data[cond_b].values

		# Perform Wilcoxon signed-rank test
		stat, p_val = wilcoxon(values_a, values_b)
		results_table.append({
			'ROI': roi,
			'Statistic': stat,
			'P_value_raw': p_val,
			'P_value_adj': np.minimum(p_val * n_obs, 1.0), # Bonferroni correction
			'N_observations': n_obs,
			'Comparison': f"{cond_a} vs {cond_b}"
		})

# Create the results dataframe
results_df = pd.DataFrame(results_table)
print(results_df.head())

# Save dataframe
results_df.to_csv(
	c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-paired_betas.csv",
	index=False
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04b. Plotting the paired results
# Another plot against zero!
# TODO: Statistics Plot per RUN: violin plots with ind. run dots connected
# TODO: Statistics Plot per VOXEL: matrix
# colors = ['#EEDF5A', '#A8C6FF']

# Figure layout
# n_rows = int(len(c.rois) / plotConf["cols"])
# fig, axes = plt.subplots(
# 	n_rows,
# 	plotConf["cols"],
# 	figsize=plotConf["figsize"],
# 	sharey=True
# )
# axes = axes.flatten()

# # Iterate through the rois
# for i, roi in enumerate(c.rois):
# 	ax = axes[i]
	
# 	# Filter data for current ROI
# 	roi_df = df[df['ROI'] == roi]
	
# 	# Create violin shape
# 	sns.violinplot(
# 		data=roi_df, 
# 		x='Condition', 
# 		y='Value',
# 		hue="Condition",
# 		palette="Set2", 
# 		order=conditions,
# 		ax=ax,
# 		inner=None,
# 		linewidth=1.5,
# 		cut=0
# 	)
	
# 	# Plot individual points and lines between them
# 	roi_pivot = roi_df.pivot(
# 		index='Idx',
# 		columns='Condition',
# 		values='Value'
# 	)
# 	rng = np.random.default_rng(42)

# 	jitters = {
# 		run: rng.uniform(-0.08, 0.08)
# 		for run in roi_pivot.index
# 	}

# 	for run in roi_pivot.index:

# 		val1 = roi_pivot.loc[run, conditions[0]]
# 		val2 = roi_pivot.loc[run, conditions[1]]

# 		j = jitters[run]

# 		x1 = 0 + j
# 		x2 = 1 + j

# 		ax.plot(
# 			[x1, x2],
# 			[val1, val2],
# 			color='gray',
# 			alpha=0.3,
# 			lw=0.8,
# 			zorder=1
# 		)

# 		ax.scatter(
# 			[x1, x2],
# 			[val1, val2],
# 			color='black',
# 			s=30,
# 			alpha=0.7,
# 			zorder=10
# 		)
	
# 	# Set labels and title
# 	ax.set_title(
# 		f'{roi}',
# 		fontsize=plotConf["subplot_fontsize"],
# 		fontweight='bold'
# 	)
# 	ax.set_xlabel("")
# 	ax.set_ylabel("")
	
# fig.suptitle(
# 	f"Averaged across voxels: {c.average_voxels}",
# 	fontsize=plotConf["fig_fontsize"],
# 	fontweight="bold"
# )
# fig.supxlabel('Condition', fontsize=plotConf["fig_fontsize"])
# fig.supylabel('Beta Estimate', fontsize=plotConf["fig_fontsize"])
# plt.tight_layout()

# if c.save_fig:
# 	fig_name = f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}-one_sample.png"
# 	fig_path = c.out_2nd / fig_name
# 	plt.savefig(fig_path, dpi=plotConf["dpi"], bbox_inches="tight")

# if c.show_fig:
# 	plt.show()
# else:
# 	plt.close(fig)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 05. Inferential statistics
# Compare betas to the residuals
# RQ: Is there a significant difference between residuals and betas?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO!