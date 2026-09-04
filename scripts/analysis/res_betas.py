#! /usr/bin/env python
# Time-stamp: <04-09-2026 m.utrosa@bcbl.eu>
"""
Extract ROI arrays from beta images and plot a single violin 
plot per ROI, where each beta is an average from all runs for that
subject (n = n_voxels) OR  across voxels per run (n = n_runs).

Before running this script ensure that you have resampled the atlas
correctly to the resolution of the functional images.
"""
# Import python packages
import pandas as pd
import numpy as np
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

# Import custom-made functions
import config as c
from config import plotConf, apply_figure_style
from utils import extract_roi_array, plot_violins_average
apply_figure_style()
# TODO: reduce the iterations => this code is slow because we're iterating
# twice through all beta images in the same way to do different things; join
# parts 02 and 01 ;)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 00. Check that inputs are defined correctly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if c.average_voxels and c.average_runs:
    raise ValueError(
        "Statistical tests cannot be performed when averaging across BOTH runs and voxels. "
        "Please set either average_voxels=False or average_runs=False."
    )

if not c.average_voxels and not c.average_runs:
    raise ValueError(
        "To perform statistical tests we need one-dimensional arrays, "
        "which means that the selected values per beta image "
        "have to be averaged EITHER across runs or voxels."
    )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Extract beta values per voxel from each ROI.
# ROI1: {C1: [[v1, v2, v3, ...],[[v1, v2, v3, ...]]], 
#        C2: [[[v1, v2, v3, ...]],[[v1, v2, v3, ...]]]} 
# Each roi is a dictionary of length n_cond
# Each cond is a list with n_run arrays
# Each run array has n_voxel values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize a dictionary to save extracted values
extracted_beta = {name: {c: [] for c in c.conditions }  for name in c.rois.keys()}
roi_names = list(extracted_beta.keys())
beta_affine = None
for sesID in c.sesIDs:                
	for acqID in c.acqIDs:

		# Inspect the SPM.xX.name to see which beta images correspond to
		# which conditions of the SPM design matrix.  Assuming numerical
		# naming of beta images: beta_space-T1wFOV_0004.nii.
		for b in range(1, len(c.conditions) + 1):
			
			# Current condition
			cond = c.conditions[b - 1]

			# Construct the path
			beta_name = f"{c.beta_filename}_{b:04d}.nii"
			beta_fold = c.dataPath / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
			beta_path = beta_fold / beta_name

			# Extract the subcortical arrays		
			mask_subcor, _, beta_subcor_affine = extract_roi_array(
				c.subID,
				sesID,
				acqID,
				c.atlas_subcor_path,
				c.space,
				beta_path,
				c.rois_subcortical,
				c.out_1st,
				verbose=False,
				save=c.save_roi,
				average_voxels=False # Keeping this false for consistency
			)
			
			# Extract the cortical arrays
			mask_cor, _, beta_cor_affine = extract_roi_array(
				c.subID,
				sesID,
				acqID,
				c.atlas_cor_path,
				c.space,
				beta_path,
				c.rois_cortical,
				c.out_1st,
				verbose=False,
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
	for cond in c.conditions:
		print(f"{roi_name} - {cond}: {np.shape(extracted_beta[roi_name][cond])}")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Transform the extracted beta values.
# Average across runs or voxels: (n_runs, n_voxels)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
selected_betas = {name: {c: [] for c in c.conditions }  for name in c.rois.keys()}

# Iterate through the ROIs
for roi_name in c.rois.keys():

	roi_dict = extracted_beta[roi_name]

	# Iterate through timing deviancy conditions	
	for condition in c.conditions:

		# Get a list of arrays (n_runs length)	
		array_list = roi_dict[condition]
		
		# Optionally, removing empty arrays	(zero values)
		if c.remove_empty:
			valid_arrays = [arr for arr in array_list if arr.size > 0]
			if len(valid_arrays) != len(array_list):
				message = f"For ROI {roi_name}, removing {len(array_list) - len(valid_arrays)} empty array(s) before averaging."
				warnings.warn(message)
		else:
			valid_arrays = array_list

		if len(valid_arrays) > 0:

			# Collapse voxels: get a mean contrast value across voxels
			if c.average_voxels:
				averaged_array = np.mean(array_list, axis=1)
				selected_betas[roi_name][condition].append(averaged_array)

			# Collapse runs: get a mean contrast value across runs
			elif c.average_runs:
				averaged_array = np.mean(array_list, axis=0)
				selected_betas[roi_name][condition].append(averaged_array)
			
			# Optionally save all averaged betas to disk
			# Why here averaged and the contrasts summed?
			if c.save_averaged:
				summed_filename = f"roi-{roi_name}_sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_job-{c.jobName}_cond-{condition}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_{c.beta_filename}.nii.gz"
				summed_path = c.out_1st / summed_filename
				
				# Save the summed array
				if beta_affine is not None:
					nib.save(nib.Nifti1Image(averaged_array, beta_affine), summed_path)

# Print update on the structure of array
print(
    "\n--- TRANSFORMED EXTRACTED DATA ---",
    f"\nAveraged across runs: {c.average_runs}",
    f"\nAveraged across voxels: {c.average_voxels}\n")
for roi_name in roi_names:
	for cond in c.conditions:
	    if np.isscalar(selected_betas[roi_name]):
	        print(f"{roi_name} - {cond}: {float(selected_betas[roi_name][cond]):.4f}")
	    else:
	        print(f"{roi_name} - {cond}: {np.shape(np.array(selected_betas[roi_name][cond]))}")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03. Descriptive plotting
# Plot beta values per ROI (subplot) and condition (x axis) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_violins_average(
	selected_betas,
	c.subID,
	c.sessions,
	c.blocks,
	c.plot_rois,
	plotConf["cols"],
	c.out_1st,
	c.space,
	scale=False,
	save=c.save_fig,
	average_runs=c.average_runs,    # only for the filename
	average_voxels=c.average_voxels # only for the filename
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04. Statistics #TODO: CHECK CORRECTNESS -- smth IS wrong : Bonferroni?!
# Compare conditions within an ROI (with averaged or not voxels).
# RQ: Is there a significant difference between conditions in the ROI?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if c.average_runs:
    data_rows = []
    
    # 1. Reshape data to keep Run x Voxel pairing
    for roi_name, cond_dict in extracted_beta.items():
        
        # Get the number of runs (should be consistent across conditions)
        # Assuming cond_dict[cond] is a list of arrays (one array per run)
        n_runs = len(cond_dict[c.conditions[0]])
        
        # Get the number of voxels (assumed constant across runs)
        # Extract from the first run of the first condition
        n_voxels = len(cond_dict[c.conditions[0]][0])
        
        for run_idx in range(n_runs):
            for voxel_idx in range(n_voxels):
                # Extract value for Condition 1
                val_cond1 = cond_dict[c.conditions[0]][run_idx][voxel_idx]
                
                # Extract value for Condition 2
                val_cond2 = cond_dict[c.conditions[1]][run_idx][voxel_idx]
                
                data_rows.append({
                    'ROI':       roi_name,
                    'Voxel_ID':  voxel_idx,    # 0-indexed voxel ID
                    'Run_ID':    run_idx,      # 0-indexed run ID
                    'Condition': c.conditions[0],
                    'Value':     val_cond1
                })
                
                data_rows.append({
                    'ROI':       roi_name,
                    'Voxel_ID':  voxel_idx,
                    'Run_ID':    run_idx,
                    'Condition': c.conditions[1],
                    'Value':     val_cond2
                })

    df = pd.DataFrame(data_rows)

    # 2. Pivot to create paired columns
    # Index: ROI + Voxel + Run. Columns: Conditions.
    pivot_df = df.pivot_table(
        index=['ROI', 'Voxel_ID', 'Run_ID'],
        columns='Condition',
        values='Value'
    )

    # 3. Perform Wilcoxon test for each voxel in each ROI
    results_table = []
    
    # Group by ROI and Voxel
    grouped = pivot_df.groupby(level=['ROI', 'Voxel_ID'])
    
    for (roi, voxel_id), group in grouped:

        # Check if both conditions exist for this group
        if c.conditions[0] in group.columns and c.conditions[1] in group.columns:
            cond1_vals = group[c.conditions[0]].values
            cond2_vals = group[c.conditions[1]].values
            
            try:
                stat, p_val = wilcoxon(cond1_vals, cond2_vals)
                results_table.append({
                    'ROI': roi,
                    'Voxel_ID': voxel_id,
                    'Statistic': stat,
                    'P-Value': p_val,
                    'N': len(cond1_vals)
                })
            except Exception as e:
                # Handle cases with constant values or insufficient data
                results_table.append({
                    'ROI': roi,
                    'Voxel_ID': voxel_id,
                    'Statistic': np.nan,
                    'P-Value': np.nan,
                    'N': len(cond1_vals)
                })
        else:
            results_table.append({
                'ROI': roi,
                'Voxel_ID': voxel_id,
                'Statistic': np.nan,
                'P-Value': np.nan,
                'N': len(cond1_vals)
            })

    results_df = pd.DataFrame(results_table)
    
    # Save dataframe
    results_df.to_csv(
        c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-paired_betas.csv",
        index=False
    )

elif c.average_voxels:
	
	# Reshape the data for plotting	
	data_rows = []
	for roi_name, cond_dict in selected_betas.items():
		for cond in c.conditions:
			values = cond_dict[cond][0]
			
			for idx, val in enumerate(values):
				data_rows.append({
					'ROI': roi_name,
					'Condition': cond,
					'Value': val,
					'Run': idx + 1,
					'N': len(values)
				})
	df = pd.DataFrame(data_rows)

	# Perform the statistical test
	results_table = []

	# Pivot the data by RUN and ROI
	pivot_df = df.pivot_table(index=['ROI', 'Run'], columns='Condition', values='Value')
	
	# Perform Wilcoxon signed-rank test
	rois = df['ROI'].unique()
	for roi in rois:
		roi_data = pivot_df.loc[roi]
		stat, p_val = wilcoxon(roi_data[c.conditions[0]], roi_data[c.conditions[1]])
		results_table.append({
			'ROI': roi,
			'Statistic': stat,
			'P-Value': p_val,
			 'N': len(values) # number of observations compared in the t-test
		})

	results_df = pd.DataFrame(results_table)

	# Save dataframe
	results_df.to_csv(
    	c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-paired_betas.csv",
    	index=False
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 05a. Statistics Plot per RUN
# TODO: Statistics Plot per VOXEL matrix
if c.average_voxels:
	sns.set_context("paper", font_scale=1.3)
	sns.set_style("white")
	plt.rcParams['font.family'] = 'sans-serif'
	plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
	plt.rcParams['axes.linewidth'] = 1.2

	colors = ['#EEDF5A', '#A8C6FF']
	hue_order = c.conditions

	n_cols = plotConf["cols"]
	n_rows = int(len(rois) / n_cols)
	fig, axes = plt.subplots(
		n_cols,
		n_rows,
		figsize=plotConf["figsize"],
		sharey=True
	)
	axes = axes.flatten()

	for i, roi in enumerate(rois):
		ax = axes[i]
		
		# Filter data for current ROI
		roi_df = df[df['ROI'] == roi]
		
		# Create violin shape
		sns.violinplot(
			data=roi_df, 
			x='Condition', 
			y='Value',
			hue="Condition",
			palette=colors, 
			order=hue_order,
			ax=ax,
			inner=None,
			linewidth=1.5,
			cut=0
		)
		
		# Plot individual points and lines between them
		roi_pivot = roi_df.pivot(
			index='Run',
			columns='Condition',
			values='Value'
		)
		rng = np.random.default_rng(42)

		jitters = {
			run: rng.uniform(-0.08, 0.08)
			for run in roi_pivot.index
		}

		for run in roi_pivot.index:

			val1 = roi_pivot.loc[run, c.conditions[0]]
			val2 = roi_pivot.loc[run, c.conditions[1]]

			j = jitters[run]

			x1 = 0 + j
			x2 = 1 + j

			ax.plot(
				[x1, x2],
				[val1, val2],
				color='gray',
				alpha=0.3,
				lw=0.8,
				zorder=1
			)

			ax.scatter(
				[x1, x2],
				[val1, val2],
				color='black',
				s=30,
				alpha=0.7,
				zorder=10
			)
		
		# Set labels and title
		ax.set_title(
			f'{roi}',
			fontsize=plotConf["subplot_fontsize"],
			fontweight='bold'
		)
		ax.set_xlabel("")
		ax.set_ylabel("")

		# Add p-value text
		p_val = results_df[results_df['ROI'] == roi]['P-Value'].values[0]
		sig_marker = ""
		if not np.isnan(p_val):
				if p_val < 0.001:
					sig_marker = "**"
				elif p_val < 0.05:
					sig_marker = "*"
				else:
					sig_marker = "ns"
				
				ax.text(0.5, 0.95, f'p = {p_val:.3f}\n{sig_marker}', 
						transform=ax.transAxes, 
						ha='center', va='top', 
						fontsize=plotConf["subplot_fontsize"]
						)
		
	fig.suptitle(
		f"Averaged across voxels: {c.average_voxels}",
		fontsize=plotConf["fig_fontsize"],
		fontweight="bold"
	)
	fig.supxlabel('Condition', fontsize=plotConf["fig_fontsize"])
	fig.supylabel('Beta Estimate', fontsize=plotConf["fig_fontsize"])
	plt.tight_layout()

	if c.save_fig:
		fig_name = f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_averageVox-{c.average_voxels}-one_sample.png"
		fig_path = c.out_2nd / fig_name
		plt.savefig(fig_path, dpi=plotConf["dpi"], bbox_inches="tight")
		plt.close(fig)
	
	if c.show_fig:
		plt.show()