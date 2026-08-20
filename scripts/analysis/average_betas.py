#! /usr/bin/env python
# Time-stamp: <24-06-2026 m.utrosa@bcbl.eu>
"""
Extract ROI arrays from beta images and plot a single violin 
plot per ROI, where each beta is an average from all runs for that
subject (n = n_voxels) OR  across voxels per run (n = n_runs).

Before running this script ensure that you have resampled the atlas
correctly to the resolution of the functional images
"""
# Which job is submitted?
# J1: whenwhat (timDev vs freqDev)
# J2: when11where (abs timDev vs freqDev)
jobName = "whenwhat"
denoising = True # NORDIC True or False

# Conditions have to be in the order of beta images
# Please see names in the SPM design matrix.
# conditions = [4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125]
conditions = ["timDev", "freqDev"]

sesIDs = [2, 3, 4, 5, 6, 7] # 2, 3, 4, 5, 6, 7
sessions = 234567 # appears in the plots' filename

acqIDs = ["BLOCK1", "BLOCK2", "BLOCK3", "BLOCK4"] # "BLOCK2", "BLOCK3", "BLOCK4"
blocks = "1234" # appears in the plots' filename

# Summation and plotting preferences
save_fig       = True  # applies to figures with statistical results 
save_roi       = False # Applies to extracted ROI arrays
save_average   = False # Applies to the summed beta arrays
average_voxels = True
remove_empty   = False # Remove or not empty arrays (e.g.: If we do not average across voxels, 
					   # do we, when averaging across runs, include voxels that have zero 
					   # beta values or not?)
# ------------------------------------------------------------------------------------------------------------------
# BELOW: DO NOT MODIFY
# ------------------------------------------------------------------------------------------------------------------
# 00. Prerequisites -------------------------------------------------
# Import python packages
import bids
import warnings
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon   # rank test

# Import custom-made functions
import grabber
import roisExtVis as rem

# 01. Parameters ----------------------------------------------------
# Experiment info
subID  = 5
anatID = 2
space  = "T1w"
plot_rois = ["IC-L", "IC-R", "MGB-L", "MGB-R", "aHG-L", "aHG-R"] 

# Project directories
homePath   = Path("/home/mutrosa/mutrosa/Documents/devLoc")
resPath    = homePath / "results"
dataPath   = resPath / jobName / f"NORDIC-{denoising}" / "1stLevel"
out_dir    = resPath / jobName / f"NORDIC-{denoising}" / "1stLevel" / "visualization"
# out_dir    = resPath / "visualization" / f"NORDIC-{denoising}" / jobName
out_dir.mkdir(parents=True, exist_ok=True)

# Get Sitek's subcortical atlas
atlas_subcor_name = f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
atlas_subcor_path = homePath / "templates" / atlas_subcor_name

# Get FreeSurfer parcellation atlas
# Desikan-Killiany Atlas: aparc+aseg.mgz
# Destrieux Atlas: aparc.a2009s+aseg.mgz (more areas)
# Tip: use fresurfers' mri_convert to convert to nii ;)
atlas_cor_name  = f"aparc.a2009s+aseg_NORDIC-{denoising}_space-{space}.nii.gz"
atlasPath       = homePath / "templates"
atlas_cor_path  =  atlasPath / atlas_cor_name

# ROIs --------------------------------------------------------------
# Size represents the volume (mm3) of subcortical structures based on
# in-vivo functional clusters (see Table 1, Sitek et al., 2019).
# Label is identified from plotting unique atlas values in freeview.
rois_subcortical = {
	# 'CN-L'  : {'size': 11, 'label': 1},
	# 'CN-R'  : {'size': 11, 'label': 2},
	# 'SOC-L' : {'size': 29, 'label': 3},
	# 'SOC-R' : {'size': 29,  'label': 4},
	'IC-L'  : {'size': 146, 'label': 5},
	'IC-R'  : {'size': 146, 'label': 6},
	'MGB-L' : {'size': 152, 'label': 7},
	'MGB-R' : {'size': 152, 'label': 8}
	}

# Fresurfer cortical areas legend
# 11133 aHG-L  G_temp_sup-G_T_transv  Anterior transverse temporal gyrus (~A1)
# 12133 aHG-R 
# 11136 PT-L   G_temp_sup-Plan_tempo  Planum temporale of the superior temporal gyrus (~A2)
# 12136 PT-R
rois_cortical = {
	'aHG-L' : {'label': 11133},
	'aHG-R' : {'label': 12133}
	}

# All rois combined
rois = rois_subcortical | rois_cortical

# 02. Extract beta values from ROIs ----------------------------------------
all_betas = {name: {c: [] for c in conditions }  for name in rois.keys()}
beta_affine = None

for sesID in sesIDs:                  
	for acqID in acqIDs:

		# Inspect the SPM.xX.name to see which beta images correspond to
		# which conditions of the SPM design matrix.  Assuming numerical
		# naming of beta images: beta_space-T1wFOV_0004.nii.
		for b in range(1, len(conditions) + 1): 
			
			# Current condition
			cond = conditions[b - 1]

			# Construct the path
			name_beta = f"beta_space-T1wFOV_{b:04d}.nii" # Have to be resampled prior to visualization!
			beta_fold = dataPath / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
			beta_path = beta_fold / name_beta

			# Extract the subcortical arrays		
			mask_subcor, _, beta_subcor_affine = rem.extract_roi_array(
				subID,
				sesID,
				acqID,
				atlas_subcor_path,
				space,
				beta_path,
				rois_subcortical,
				out_dir,
				verbose=False,
				save=save_roi,
				average_voxels=average_voxels
			)
			
			# Extract the cortical arrays
			mask_cor, _, beta_cor_affine = rem.extract_roi_array(
				subID,
				sesID,
				acqID,
				atlas_cor_path,
				space,
				beta_path,
				rois_cortical,
				out_dir,
				verbose=False,
				save=save_roi,
				average_voxels=average_voxels
			)
		
			# Accumulate subcortical arrays for summation
			for name in rois_subcortical.keys():
				all_betas[name][cond].append(mask_subcor[name])

			# Accumulate cortical arrays for summation
			for name in rois_cortical.keys():
				all_betas[name][cond].append(mask_cor[name])

# Assign beta affine
if beta_cor_affine.all() == beta_subcor_affine.all():
	beta_affine = beta_subcor_affine

# 03. Average -------------------------------------------------------
if average_voxels == False:
	summed_betas = {name: {c: [] for c in conditions }  for name in rois.keys()}

	# Iterate through the ROIs
	for roi_name in rois.keys():

		roi_dict = all_betas[roi_name]

		# Iterate through timing deviancy conditions	
		for condition in conditions:

			# Get a list of arrays (n_runs length)	
			array_list = roi_dict[condition]

			# Optionally, removing empty arrays	(zero values)
			if remove_empty:
				valid_arrays = [arr for arr in array_list if arr.size > 0]
				if len(valid_arrays) != len(array_list):
					message = f"For ROI {roi_name}, removing {len(array_list) - len(valid_arrays)} empty array(s) before averaging."
					warnings.warn(message)
			else:
				valid_arrays = array_list

			if len(valid_arrays) > 0:

				# Assuming that every array (per ROI) has the same length (must be True if ROI anatomically defined)
				# Note: len(valid_arrays) == len(sesIDs) * len(acqIDs) unless empty arrays are removed
				averaged_array = np.mean(array_list, axis=0)
				summed_betas[roi_name][condition].append(averaged_array)

				# Optionally save all averaged betas to disk
				if save_average:
					summed_filename = f"sub-{subID:02d}_roi-{roi_name}_space-T1wFOV_cond-{condition}_type-average.nii.gz"
					summed_path = out_dir / summed_filename
					
					# Save the summed array
					if beta_affine is not None:
						nib.save(nib.Nifti1Image(averaged_array, beta_affine), summed_path)

# 04. Plot ----------------------------------------------------------
# Plot beta values per ROI (subplot) and condition (x axis) 
# NOT averaged across ROI voxels: one value per voxel (average_voxels=False)
# Averaged across runs (summed_betas)
# Each violin has n_voxels dots
if average_voxels == False:
	rem.plot_violins_average(
		summed_betas,
		subID,
		sessions,
		blocks,
		plot_rois,
		2,
		out_dir,
		space,
		scale=False,
		save=True,
		average_runs=True,
		average_voxels=average_voxels)

# Plot beta values per ROI (subplot) and condition (x axis) 
# Averaged across ROI voxels: one value per ROI (average_voxels=True)
# NOT averaged across runs (all_betas)
# Each violin has n_runs dots
if average_voxels == True:
	rem.plot_violins_average(
		all_betas, # The difference is here (in respect to the previous violin plot figure)
		subID,
		sessions,
		blocks,
		plot_rois,
		2,
		out_dir,
		space,
		scale=False,
		save=True,
		average_runs=False,
		average_voxels=average_voxels)
	
# 05. Statistics ----------------------------------------------------------
# Compare conditions within an ROI: Is there a significant difference
# between conditions in the ROI? Comparing voxels within an ROI.
if average_voxels == False: #TODO 	CHECK BELOW !!
    data_rows = []
    
    # 1. Reshape data to keep Run x Voxel pairing
    for roi_name, cond_dict in summed_betas.items():
        # Get the number of runs (should be consistent across conditions)
        # Assuming cond_dict[cond] is a list of arrays (one array per run)
        n_runs = len(cond_dict[conditions[0]])
        
        # Get the number of voxels (assumed constant across runs)
        # Extract from the first run of the first condition
        n_voxels = len(cond_dict[conditions[0]][0])
        
        for run_idx in range(n_runs):
            for voxel_idx in range(n_voxels):
                # Extract value for Condition 1
                val_cond1 = cond_dict[conditions[0]][run_idx][voxel_idx]
                
                # Extract value for Condition 2
                val_cond2 = cond_dict[conditions[1]][run_idx][voxel_idx]
                
                data_rows.append({
                    'ROI': roi_name,
                    'Voxel_ID': voxel_idx,  # 0-indexed voxel ID
                    'Run_ID': run_idx,      # 0-indexed run ID
                    'Condition': conditions[0],
                    'Value': val_cond1
                })
                
                data_rows.append({
                    'ROI': roi_name,
                    'Voxel_ID': voxel_idx,
                    'Run_ID': run_idx,
                    'Condition': conditions[1],
                    'Value': val_cond2
                })

    df = pd.DataFrame(data_rows)

    # 2. Pivot to create paired columns
    # Index: ROI + Voxel + Run. Columns: Conditions.
    pivot_df = df.pivot_table(
        index=['ROI', 'Voxel_ID', 'Run_ID'],
        columns='Condition',
        values='Value'
    )

    # 3. Perform Wilcoxon test for each Voxel in each ROI
    results_table = []
    
    # Group by ROI and Voxel
    grouped = pivot_df.groupby(level=['ROI', 'Voxel_ID'])
    
    for (roi, voxel_id), group in grouped:
        # Check if both conditions exist for this group
        if conditions[0] in group.columns and conditions[1] in group.columns:
            cond1_vals = group[conditions[0]].values
            cond2_vals = group[conditions[1]].values
            
            try:
                stat, p_val = wilcoxon(cond1_vals, cond2_vals)
                results_table.append({
                    'ROI': roi,
                    'Voxel_ID': voxel_id,
                    'Statistic': stat,
                    'P-Value': p_val
                })
            except Exception as e:
                # Handle cases with constant values or insufficient data
                results_table.append({
                    'ROI': roi,
                    'Voxel_ID': voxel_id,
                    'Statistic': np.nan,
                    'P-Value': np.nan
                })
        else:
            results_table.append({
                'ROI': roi,
                'Voxel_ID': voxel_id,
                'Statistic': np.nan,
                'P-Value': np.nan
            })

    results_df = pd.DataFrame(results_table)

elif average_voxels == True:
	
	# Reshape the data for plotting	
	data_rows = []
	for roi_name, cond_dict in all_betas.items():
		conditions = list(cond_dict.keys())
		for cond in conditions:
			values = cond_dict[cond]
			
			for idx, val in enumerate(values):
				data_rows.append({
					'ROI': roi_name,
					'Condition': cond,
					'Value': val,
					'Run': idx + 1
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
		stat, p_val = wilcoxon(roi_data[conditions[0]], roi_data[conditions[1]])
		results_table.append({
			'ROI': roi,
			'Statistic': stat,
			'P-Value': p_val
		})

	results_df = pd.DataFrame(results_table)

# 06. Statistics Plot per RUN-----------------------------------------------------
if average_voxels == True:
	sns.set_context("paper", font_scale=1.3)
	sns.set_style("white")
	plt.rcParams['font.family'] = 'sans-serif'
	plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
	plt.rcParams['axes.linewidth'] = 1.2

	colors = ['#EEDF5A', '#A8C6FF']
	hue_order = conditions

	n_cols = 3
	n_rows = int(len(rois) / n_cols)
	fig, axes = plt.subplots(n_cols, n_rows, figsize=(12, 14), sharey=True)
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
			linewidth=1.5
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

			val1 = roi_pivot.loc[run, conditions[0]]
			val2 = roi_pivot.loc[run, conditions[1]]

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
		ax.set_title(f'{roi}', fontsize=11, fontweight='bold')
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
						fontsize=12
						)
		
	fig.suptitle("Averaged across voxels", fontsize=12, fontweight="bold")
	fig.supxlabel('Condition', fontsize=14)
	fig.supylabel('Beta Estimate', fontsize=14)
	plt.tight_layout()

	if save_fig:
		fig_name = f"sub-{subID:02d}_ses-{sessions}_block-{blocks}_space-{space}_job-{jobName}_averageVox-{average_voxels}-one_sample.png"
		fig_path = out_dir / fig_name
		plt.savefig(fig_path, dpi=300, bbox_inches="tight")
		plt.close(fig)
	plt.show()