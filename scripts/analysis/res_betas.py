#! /usr/bin/env python
# Time-stamp: <17-09-2026 m.utrosa@bcbl.eu>
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
from utils import extract_roi_array, plot_violins_zero_betas
from utils import plot_violins_betas_paired
# IMPORT ALL FIGURE FUNCTIONS HERE AS SEPARATE FUNCTIONS

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
# 03a. Inferential statistics: Beta distributions compared against zero
# Compare beta distributions within an ROI (either across runs or voxels).
# RQ: Is there a sig. difference of the effect of the regressor from 0?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("\nPerforming two-sided Wilcoxon rank tests against 0 for each beta.")
results_one_sample = []
for roi in roi_names:
	for cond in conditions:

		data = np.array(selected_betas[roi][cond])

		# Check dimensions of the data
		if data.ndim == 0:
			raise ValueError(
				f"\nFor ROI {roi} and regressor {cond} "
				"the data is a scalar (N=1). "
				"Cannot perform statistics. Check averaging settings.")
		elif data.ndim > 1:
			data = data.flatten() # Wilcoxon test requires one-dimensional input

		# One-sample Wilcoxon Signed-Rank Test (non-parametric)
		# Tests if the median of the distribution is different from 0
		res_wilcox = wilcoxon(
			x = data,                   # Must be 1D
			y = None,                   # None implies one-sample test against 0
			zero_method = "pratt",      # Includes zero-differences in the ranking process
			correction  = False,        # Default
			alternative = "two-sided",
			method      = "auto"        # Default
		)

		statistic = res_wilcox.statistic
		p_value   = res_wilcox.pvalue

		# Calculate metrics (and control for division by zero)
		mean_val = np.mean(data)
		std_val  = np.std(data, ddof=c.ddof)
		cohen_d  = mean_val / std_val if std_val > 0 else np.nan

		# Append to list
		results_one_sample.append({
			"ROI"            : roi,
			"regressor"      : cond,
			"N"              : len(data),
			"mean_beta"      : mean_val,
			"std_beta"       : std_val,
			"stat"           : statistic,
			"p_value"        : p_value,
			"cohen_d"        : cohen_d
		})

# Create and display DataFrame
df_results_against0 = pd.DataFrame(results_one_sample)
print(df_results_against0.head())

# Calculate the number of tests for Bonferroni correction
# Correction has to be done on the data that is tested multiple times # CHECK
n_tests = len(conditions)
print(f"\nThe Bonferroni correction is applied for {len(roi_names)} tests.")

# Adjust p-values with Bonferroni and reduce the adjusted p-values that exceed 1 to 1
df_results_against0['p_value_bonferroni'] = df_results_against0['p_value'] * n_tests
df_results_against0['p_value_bonferroni'] = df_results_against0['p_value_bonferroni'].clip(upper=1.0)

# Save dataframe
df_results_against0.to_csv(
    c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-against0-betas.csv",
    index=False
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03b. Inferential statistics: Beta pairs per ROI
# Compare beta distributions within an ROI (either across runs or voxels).
# RQ: Is there a sig. difference between the effects of regressor pairs?
# Filter for betas that show a significant effect from zero ? # CHECK
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
				'nVox_nRun': selected_betas[roi_name][cond][0].shape,
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

# Pivot the data by index (run or voxel) and ROI
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

		# Number of tests and observations
		n_tests = len(condition_pairs)
		n_obs = len(roi_data)

		# Extract the arrays
		values_a = valid_pair_data[cond_a].values
		values_b = valid_pair_data[cond_b].values

		# Perform Wilcoxon signed-rank test
		stat, p_val = wilcoxon(values_a, values_b)
		results_table.append({
			'ROI': roi,
			'Statistic': stat,
			'P_value_raw': p_val,
			'P_value_adj': np.minimum(p_val * n_tests, 1.0), # Bonferroni correction
			'N_observations': n_obs,
			'Comparison': [cond_a, cond_b]
		})

# Create the results dataframe
df_results_paired = pd.DataFrame(results_table)
print(df_results_paired.head())

# Save dataframe
df_results_paired.to_csv(
	c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-paired_betas.csv",
	index=False
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04a. Figure 1: Plotting the against-zero results
# AIM: see which ROIs distinguish between which betas (paired)
# RQ: Is the regressor effect sig. different from 0?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_violins_zero_betas(
		selected_betas,
		df_results_against0,
		c.subID,
		c.plot_rois,
		c.plotConf,
		c.out_2nd,
		c.space, # only for the filename
		c.jobName,
		save=c.save_fig,
		show=c.show_fig,
		average_runs=c.average_runs,    # only for the filename
		average_voxels=c.average_voxels # only for the filename
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04b. Figure 2 (version runs): Plotting the paired results
# AIM: see which ROIs (average voxels) distinguish between which betas (paired)
# RQ: Which regressor pairs do ROIs distinguish significantly?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if c.average_voxels:
	
	plot_violins_betas_paired(
		selected_betas,
		df_results_paired,
		c.subID,
		c.out_2nd,
		c.space,   # only for the filename
		c.jobName, # only for the filename
		c.plotConf,
		save=c.save_fig,
		show=c.show_fig,
		average_runs=c.average_runs,    # only for the filename
		average_voxels=c.average_voxels # only for the filename
	)
	# # Filter for sig. results
	# p_value = "P_value_raw"
	# filtered_df = df_results_paired[df_results_paired[p_value] < 0.05]
	# filtered_df.reset_index(drop=True, inplace=True)

	# if filtered_df.empty:
	# 	print("No ROI distinguishes between a regressor pair significantly.")
	# else:

	# 	# Figure layout
	# 	fig, axes = plt.subplots(
	# 		int(np.ceil(len(filtered_df) / plotConf["cols"])),
	# 		plotConf["cols"],
	# 		figsize=plotConf["figsize"],
	# 		sharey=True
	# 	)
	# 	axes = axes.flatten()

	# 	# Colors
	# 	cb_palette = ["#FFC20A", "#0C7BDC"]

	# 	# Iterate through the filtered rois
	# 	for row in filtered_df.itertuples():
	# 		ax = axes[row.Index]

	# 		# Get conditions
	# 		cond1 = row.Comparison[0]
	# 		cond2 = row.Comparison[1]

	# 		# Extract the sig. regressor pair data
	# 		roi_df1 = selected_betas[roi][cond1][0]
	# 		roi_df2 = selected_betas[roi][cond2][0]
	# 		n_runs  = len(roi_df2)

	# 		# Create plot df
	# 		plot_df = pd.DataFrame({
	# 			'ROI': row.ROI,
	# 			'Cond': [cond1] * n_runs + [cond2] * n_runs,
	# 			'Value': np.concatenate([roi_df1, roi_df2]),
	# 			'Run_ID': list(range(1, n_runs + 1)) * 2
	# 			})

	# 		# Ensure consistency in coloring
	# 		unique_rois = plot_df['ROI'].unique()
	# 		current_palette = {cond1: cb_palette[0], cond2: cb_palette[1]}

	# 		# A. Plot split violins
	# 		sns.violinplot(
	# 			data=plot_df, 
	# 			x='Cond', 
	# 			y='Value',
	# 			hue="Cond",
	# 			palette=current_palette, 
	# 			hue_order=[cond1, cond2],
	# 			split=False, 
	# 			inner=None,
	# 			linewidth=1.5,
	# 			alpha=0.8,
	# 			legend=False,
	# 			cut=0,
	# 			ax=ax)
			
	# 		# B. Plot individual points with controlled jitter + lines between dots
	# 		rng = np.random.default_rng(42)
	# 		jitter_strength = 0.12

	# 		# Deterministic jitter for pairing consistency
	# 		jitter_a = rng.uniform(-jitter_strength, jitter_strength, n_runs)
	# 		jitter_b = rng.uniform(-jitter_strength, jitter_strength, n_runs)
			
	# 		x_a = np.zeros(n_runs) + jitter_a
	# 		x_b = np.ones(n_runs) + jitter_b

	# 		# Scatter individual constrast estimate points
	# 		ax.scatter(x_a, roi_df1,
	# 				color='black', s=35, alpha=0.8,
	# 				edgecolor='black', linewidth=1, zorder=10)
	# 		ax.scatter(x_b, roi_df2,
	# 				color='black', s=35, alpha=0.8,
	# 				edgecolor='black', linewidth=1, zorder=10)

	# 		# Draw paired connections
	# 		for i in range(n_runs):
	# 			ax.plot([x_a[i], x_b[i]],
	# 					[roi_df1[i], roi_df2[i]],
	# 					color='gray', linewidth=1, alpha=0.6,
	# 					zorder=5, solid_capstyle='round')

	# 		# Spines
	# 		ax.spines['top'].set_visible(False)
	# 		ax.spines['right'].set_visible(False)
	# 		ax.spines['bottom'].set_linewidth(1)
	# 		ax.spines['left'].set_linewidth(1)

	# 		# Annotate significance
	# 		p_val = row.P_value_raw
	# 		is_sig = p_val < 0.05

	# 		# Get y limits
	# 		y_min, y_max = ax.get_ylim()
	# 		y_range = y_max - y_min
	# 		star_y = y_max + (y_range * 0.02)
	# 		text_y = y_max - (y_range * 0.02)
	# 		bracket_y = y_max - (y_range * 0.01)

	# 		if is_sig:
	# 			ax.plot([0, 1], [bracket_y, bracket_y], color='black', linewidth=1.4)
	# 			ax.plot([0, 0], [bracket_y, bracket_y - (y_range*0.02)], color='black', linewidth=1.4)
	# 			ax.plot([1, 1], [bracket_y, bracket_y - (y_range*0.02)], color='black', linewidth=1.4)
	# 			ax.text(0.5, star_y, '***', ha='center', va='bottom', fontsize=plotConf["fig_fontsize"] - 2, color='black')

	# 			if p_val < 0.001:
	# 				label = '***'
	# 				p_str = f"{p_val:.3f}"
	# 			elif p_val < 0.01:
	# 				label = '**'
	# 				p_str = f"{p_val:.3f}"
	# 			else:
	# 				label = '*'
	# 				p_str = f"{p_val:.3f}"
	# 			color = 'black'
	# 			font_weight = "bold"
	# 		else:
	# 			label = "ns"
	# 			p_str = ""
	# 			color = 'gray'
	# 			font_weight = "normal"
			
	# 		# Set labels and title
	# 		ax.set_title(
	# 			f'{row.ROI}',
	# 			y=1.05,
	# 			fontsize=plotConf["subplot_fontsize"],
	# 			fontweight='bold')
	# 		ax.set_xlabel("")
	# 		ax.set_ylabel("")
			
	# 	# Hide empty subplots
	# 	for j in range(len(filtered_df), len(axes)):
	# 		fig.delaxes(axes[j])

	# 	# Titles and axis lables
	# 	fig.suptitle(f"sub-{c.subID:02d}, avgVox: {c.average_voxels}, avgRun: {c.average_runs}",
	# 				 fontsize=plotConf["fig_fontsize"],
	# 				 fontweight="bold")
	# 	fig.supxlabel('Timing Deviation [msec]', fontsize=plotConf["fig_fontsize"], fontweight="bold")
	# 	fig.supylabel('Beta Estimate [β]', fontsize=plotConf["fig_fontsize"], fontweight="bold")
	# 	plt.tight_layout()

	# 	if c.save_fig:
	# 		fig_name = f"sub-{c.subID:02d}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}-paired-betas.png"
	# 		fig_path = c.out_2nd / fig_name
	# 		plt.savefig(fig_path, dpi=plotConf["dpi"], bbox_inches="tight")

	# 	if c.show_fig:
	# 		plt.show()
	# 	else:
	# 		plt.close(fig)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04c. Figure 2 (version voxels): Plotting the paired results with surface plot
# AIM: see which ROI voxels distinguish between which betas (paired)
# TODO: nilearn.plotting.plot_stat_map / plot_surf_stat_map
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: identify significant clusters
# TODO: plot significant clusters (size and location on a 3D brain)
# TODO: add a line around the actual ROI (to see how much of it was not sig.)
# Issue: how to identify the regressors well -- 22 vs 11 colors?