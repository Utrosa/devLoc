#! /usr/bin/env python
# Time-stamp: <22-06-2026 m.utrosa@bcbl.eu>
"""
2nd Level Analysis: statistical analysis on contrast img
Returns:
- effect sizes, p-values, and t-tests for each region pair
- figure showing whether an ROI pair differs in its response to a specific contrast
"""
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import python packages
import bids
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon   # rank test    
from itertools import combinations # generate all possible combos
                                   # order does not matter 

# Import custom-made functions
import grabber
import roisExtVis as rem

# 00. Prerequisites -------------------------------------------------
# *~*~*~*~ Which job is submitted? *~*~*~*~
# J1: whenwhat (timDev vs freqDev)
# J2: when11where (abs timDev vs freqDev)
jobName = "whenwhat"
denoising = True # NORDIC True or False

if jobName == "whenwhat":
    conditions = ["timDev", "freqDev"]
elif jobName == "when11where":
    conditions = [4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125]

# *~*~*~*~ Summation and plotting preferences *~*~*~*~
save_roi       = False # Applies to extracted ROI arrays
save_fig       = True
save_average   = False # Applies to the summed contrast arrays
average_voxels = True  # Do we average voxels per ROI or not? 
plot_rois      = ["IC-L", "IC-R", "MGB-L", "MGB-R", "aHG-L", "aHG-R"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BELOW DO NOT MODIFY
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *~*~*~*~ Project directories *~*~*~*~
homePath  = Path("/home/mutrosa/mutrosa/Documents/devLoc")
resPath   = homePath / "results"
dataPath  = resPath / jobName / f"NORDIC-{denoising}" / "1stLevel"
outPath   = resPath / jobName / f"NORDIC-{denoising}" / "2stLevel"
outPath.mkdir(parents=True, exist_ok=True)

# *~*~*~*~ Experiment info *~*~*~*~
subID  = 5
anatID = 2
space  = "T1w"
sesIDs = [2, 3, 4, 5, 6, 7] # 2, 3, 4, 5, 6, 7
sessions = 234567 # appears in the filenames
acqIDs = ["BLOCK1", "BLOCK2", "BLOCK3", "BLOCK4"]
blocks = "1234" # appears in the filenames

# *~*~*~*~ Atlases *~*~*~*~
# Get Sitek's subcortical atlas
atlas_subcor_name = f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
atlas_subcor_path = homePath / "templates" / atlas_subcor_name

# Get FreeSurfer's parcellation: Destrieux Atlas
atlas_cor_name  = f"aparc.a2009s+aseg_NORDIC-{denoising}_space-{space}.nii.gz"
atlasPath       = homePath / "templates"
atlas_cor_path  =  atlasPath / atlas_cor_name

# *~*~*~*~ ROIs *~*~*~*~
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
rois_cortical = {
	'aHG-L' : {'label': 11133},
	'aHG-R' : {'label': 12133}
	}

# All rois combined
rois = rois_subcortical | rois_cortical

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Average all contrast images into one and save to disk ---------
contrast_info = []
for sesID in sesIDs:
    for acqID in acqIDs:
        con_dict = {}; con_dict["sesID"] = sesID; con_dict["acqID"] = acqID
        name_contrast = "con_space-T1wFOV_0001.nii" # T1 > boldref FOV
        contrast_fold = dataPath / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
        contrast_path = contrast_fold / name_contrast
        con_dict["con_path"] = contrast_path
        contrast_info.append(con_dict)

con_img = sum([nib.load(con["con_path"]).get_fdata() for con in contrast_info])
contrast_affine = nib.load(contrast_info[0]["con_path"]).affine
if save_average:
    nib.save(
            nib.Nifti1Image(con_img, contrast_affine),
            outPath / f"con_ses-{sessions}_acq-BLOCK{blocks}_0001.nii.gz"
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Extract contrast values per voxel from each ROI ---------------
# If average_voxels == False, return for each ROI, a list of n_runs
# arrays with n_voxel values. If True, one value (array) per ROI.
# Shape: (n_runs, n_voxels)
all_cons = {name : [] for name in rois.keys()}
for contrast_dict in contrast_info:
        
        # Extract the subcortical arrays		
		mask_subcor, _, con__subcor_affine = rem.extract_roi_array(
            subID,
            contrast_dict["sesID"],
            contrast_dict["acqID"],
            atlas_subcor_path,
            space,
            contrast_dict["con_path"],
            rois_subcortical,
            outPath,
            verbose=False,
            save=save_roi,
            average_voxels=average_voxels 
		)
		
        # Extract the cortical arrays
		mask_cor, _, con_cor_affine = rem.extract_roi_array(
            subID,
            contrast_dict["sesID"],
            contrast_dict["acqID"],
            atlas_cor_path,
            space,
            contrast_dict["con_path"],
            rois_cortical,
            outPath,
            verbose=False,
            save=save_roi,
            average_voxels=average_voxels
        )
		
        # Accumulate subcortical arrays for summation
		for name in rois_subcortical.keys():
			all_cons[name].append(mask_subcor[name])

        # Accumulate cortical arrays for summation
		for name in rois_cortical.keys():
			all_cons[name].append(mask_cor[name])

# Invert structure of the dict 
# For each ROI, have an array of shape (n_voxels, n_runs)
all_cons_trans = {}
for name, array_list in all_cons.items():
    all_cons_trans[name] = np.stack(array_list, axis=1)

# Get mean contrast values per each run
mean_con = {}
for name, array_list in all_cons_trans.items():
    mean_con[name] = np.mean(array_list, axis=0) #TODO: include empty or not? 
roi_names = list(mean_con.keys())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03a. Compare contrasts per ROI against ZERO -----------------------
# Question: Does this specific ROI distinguish between freqDev and timDev?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# One sample t-tests
results_one_sample = []
for roi in roi_names:
    data = mean_con[roi]
    
    # Remove NaNs if any (optional but safe) #TODO: include empty or not? 
    data_clean = data[~np.isnan(data)]

    # One-sample Wilcoxon Signed-Rank Test (non-parametric)
    # Tests if the median of the distribution is different from 0
    # regardless of the distibution of the data.
    res_wilcox = wilcoxon(   # TODO check this function
        x = data_clean, 
        y = None,                  # None implies one-sample test against 0
        zero_method = "pratt", 
        correction  = False, 
        alternative = "two-sided",
        method = "auto"
    )
    
    statistic = res_wilcox.statistic
    p_value   = res_wilcox.pvalue
    
    # Calculate effect size
    cohen_d = np.mean(data_clean) / np.std(data_clean, ddof=1) if np.std(data_clean, ddof=1) > 0 else 0
    
    results_one_sample.append({
        "roi"       : roi,
        "n_runs"    : len(data_clean),
        "mean_contrast": np.mean(data_clean),
        "std_contrast" : np.std(data_clean, ddof=1),
        "stat"      : statistic,
        "p_value"   : p_value,
        "cohen_d"   : cohen_d
    })

# Create DataFrame and apply Bonferroni correction
df_results_against0 = pd.DataFrame(results_one_sample)
n_tests = len(df_results_against0)

# Adjust p-values
df_results_against0['p_value_bonferroni'] = df_results_against0['p_value'] * n_tests
df_results_against0['p_value_bonferroni'] = df_results_against0['p_value_bonferroni'].clip(upper=1.0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03b. Compare pairs of ROIs to see if they respond differently ------
# Is the contrast estimate larger in one ROI than the other across runs?
# Paired rank sum tests
results_paired = []

for roi_a, roi_b in combinations(roi_names, 2):
    data_a = mean_con[roi_a]
    data_b = mean_con[roi_b]
    
    # Ensure lengths match (they should if n_runs is consistent)
    if len(data_a) != len(data_b):
        continue
        
    # The Wilcoxon rank-sum test
    # The null hypothesis: the two related paired samples (rois) come from the same distribution. 
    # See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html#scipy.stats.wilcoxon
    diff_data = np.round(data_a - data_b, decimals=8)
    res = wilcoxon(
        x = diff_data,             # Differences between two sets of measurements (1D)  
        y = None,                  # Keep undefined to avoid roundoff error 
        zero_method = "pratt",     # Conservative: includes zero-differences, but drops ranks of the zeros 
                                   # DOI:10.1080/01621459.1959.10501526 and DOI:10.2307/3001968
        correction  = False,       # Default: does not apply continuinty correction
        alternative = "two-sided", # d = x 
        method = "auto",           # How to compute p-value: auto, exact, asymptotic
        axis = 0,
        nan_policy = "propagate"   # How to handle nan-values: propagate or omit 
    )
    results_paired.append({
        "roi_a"   : roi_a,
        "roi_b"   : roi_b,
        "stat"    : res.statistic, # The sum of the ranks of the differences above or below zero (for two-sided)
        "p_value" : res.pvalue,
        "cohen_d" : np.mean(diff_data) / np.std(diff_data, ddof = 1) # ddof is 1 for sample SD (Bessel's correction)    

    })

# Create a dataframe and apply Bonferroni correction
df_results_paired = pd.DataFrame(results_paired)
n_tests    = len(df_results_paired)

# Adjust p-vaues with Bonferroni correction
df_results_paired['p_value_bonferroni'] = df_results_paired['p_value'] * n_tests

# Reduce the adjusted p-values that exceed 1 to 1
df_results_paired['p_value_bonferroni'] = df_results_paired['p_value_bonferroni'].clip(upper=1.0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4a. Plotting -------------------------------------------------------
# RQ: Are the contrast values within each individual ROI 
# significantly different from zero?
# *~*~*~*~ Style settings
sns.set_context("paper", font_scale=1.3)
sns.set_style("white")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 1.2

# Define a single color for the contrast freqDev - timDev
contrast_color = "#9EDDFB" 

# *~*~*~*~ Setup grid and axes
n_rois = len(df_results_against0)
cols = 2
rows = int(np.ceil(n_rois / cols))

fig, axes = plt.subplots(
    rows,
    cols,
    figsize=(12,15),
    constrained_layout=True
)
axes = axes.flatten()

# *~*~*~*~ Sort the dataframe
df_plot = df_results_against0.copy()
df_plot = df_plot.sort_values('cohen_d', key=abs, ascending=False).reset_index(drop=True)

# *~*~*~*~ Plotting loop
for idx, (_, row) in enumerate(df_plot.iterrows()):
    ax = axes[idx]
    roi_name = row['roi']
    
    # Get raw data for this specific ROI
    data = mean_con[roi_name]
    data_clean = data[~np.isnan(data)]
    n_runs = len(data_clean)
    
    # Create a dataframe for seaborn
    plot_df = pd.DataFrame({
        'ROI': [roi_name] * n_runs,
        'Value': data_clean,
        'Run_ID': list(range(1, n_runs + 1))
    })
    
    # A. Plot violin
    sns.violinplot(
        data=plot_df, 
        x='ROI', 
        y='Value',
        ax=ax,
        color=contrast_color, 
        inner=None,          # No internal quartiles to keep it clean
        linewidth=1.2,
        width=0.5
    )
    
    # B. Plot individual points (the 24 runs)
    sns.swarmplot(
        data=plot_df,
        x='ROI',
        y='Value',
        ax=ax,
        color='black',
        size=6,
        zorder=10
    )

    # C. Overlay Mean and Error Bar (SEM)
    mean_val = np.mean(data_clean)
    sem_val = np.std(data_clean, ddof=1) / np.sqrt(n_runs)
    
    # Plot mean as a large white dot with black edge
    ax.scatter(0.03, mean_val, 
               color="#F7F3EF", edgecolor='#F7F3EF', s=35, marker='o', 
               zorder=11, linewidth=0.6, label='Mean')
    
    # Plot Error Bar (SEM)
    ax.errorbar(0.03, mean_val, yerr=sem_val, 
                fmt='none', ecolor='#FF8C00', linewidth=1.2, 
                capsize=6, capthick=2, zorder=12)
    
    # Reference line at ZERO (Crucial for one-sample test)
    ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.9, zorder=1)

    # Significance Annotation (Against Zero)
    p_val = row['p_value_bonferroni']
    is_sig = p_val < 0.05
    
    # Determine position for stars (above the max value or above the error bar)
    y_max = ax.get_ylim()[1]
    y_range = y_max - ax.get_ylim()[0]
    
    # Dynamic positioning: place stars above the data
    star_y = y_max + (y_range * 0.05) 
    text_y = y_max + (y_range * 0.08)
    
    # Adjust ylim to make room for stars if necessary
    current_ylim = ax.get_ylim()
    if star_y > current_ylim[1]:
        ax.set_ylim(current_ylim[0], star_y + (y_range * 0.05))

        # Recalculate positions after ylim change
        star_y = ax.get_ylim()[1] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05
        text_y = ax.get_ylim()[1] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02

    if is_sig:
        if p_val < 0.001:
            stars = '***'
            p_str = f"{p_val:.1e}"
        elif p_val < 0.01:
            stars = '**'
            p_str = f"{p_val:.4f}"
        else:
            stars = '*'
            p_str = f"{p_val:.4f}"
        
        ax.text(0, star_y, stars, ha='center', va='bottom', 
                fontsize=11, fontweight ="bold", color='black')
        ax.text(0, text_y, p_str, ha='center', va='bottom', 
                fontsize=12, color='black')
    else:
        ax.text(0, star_y, "ns", ha='center', va='bottom', 
                fontsize=12, color='gray', style='italic')
    
    # Title shows Mman contrast and effect size
    ax.set_title(
        f"{roi_name} | μ={mean_val:.3f}, d={row['cohen_d']:.2f}",
        fontsize=12,
        pad=8
    )
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([])

    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)

# Hide empty subplots
for j in range(n_rois, len(axes)):
    fig.delaxes(axes[j])

# Global y label
fig.text(
    0,
    0.5,
    'Contrast Estimate (freqDev - timDev)',
    va='center',
    rotation='vertical',
    fontsize=14
)

# *~*~*~*~ Optionally save
if save_fig:
    fig_name = f"sub-{subID:02d}_ses-{sessions}_block-{blocks}_space-{space}_contrast-{jobName}_one_sample.png"
    fig_path = outPath / fig_name
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

# *~*~*~*~ Show
plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Plotting -------------------------------------------------------
# RQ: 
# Plots mean contrast estimates (array of 24 values) for each ROI pair
# *~*~*~*~ Style settings
sns.set_context("paper", font_scale=1.3)
sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 1.2
cb_palette = ["#FFBB28", "#8B84D8"]

# *~*~*~*~ Setup grid and axes
cols = 3
rows = int(np.ceil(n_tests / cols))
fig, axes = plt.subplots(rows, cols, figsize=(12,14), constrained_layout=True)
axes = axes.flatten()

# *~*~*~*~ Sort the dataframe
df_plot = df_results_paired.copy()
df_plot = df_plot.sort_values('cohen_d', ascending=False).reset_index(drop=True)

# *~*~*~*~ Plotting loop
for idx, (_, row) in enumerate(df_plot.iterrows()):
    ax = axes[idx]

    # Get roi name
    roi_a_name = row['roi_a']
    roi_b_name = row['roi_b']
    
    # Get raw data
    data_a = mean_con[roi_a_name]
    data_b = mean_con[roi_b_name]
    n_runs = len(data_a)
    
    # Create a dataframe for seaborn
    plot_df = pd.DataFrame({
        'ROI': [roi_a_name] * n_runs + [roi_b_name] * n_runs,
        'Value': np.concatenate([data_a, data_b]),
        'Run_ID': list(range(1, n_runs + 1)) * 2
    })
    
    # Ensure consistency in coloring
    unique_rois = plot_df['ROI'].unique()
    current_palette = {roi_a_name: cb_palette[0], roi_b_name: cb_palette[1]}

    # A. Plot split violins
    sns.violinplot(
        data=plot_df, 
        x='ROI', 
        y='Value',
        hue='ROI', 
        ax=ax,
        palette=current_palette, 
        hue_order=[roi_a_name, roi_b_name],
        split=True, 
        inner=None,
        linewidth=1.5,
        alpha=0.8,
        legend=False
    )
    
    # B. Plot individual points with controlled jitter + correct connections
    rng = np.random.default_rng(42)
    jitter_strength = 0.12

    # deterministic jitter for pairing consistency
    jitter_a = rng.uniform(-jitter_strength, jitter_strength, n_runs)
    jitter_b = rng.uniform(-jitter_strength, jitter_strength, n_runs)

    x_a = np.zeros(n_runs) + jitter_a
    x_b = np.ones(n_runs) + jitter_b

    # Scatter individual constrast estimate points
    ax.scatter(x_a, data_a,
            color='black', s=35, alpha=0.8,
            edgecolor='black', linewidth=1, zorder=10)

    ax.scatter(x_b, data_b,
            color='black', s=35, alpha=0.8,
            edgecolor='black', linewidth=1, zorder=10)

    # Draw paired connections
    for i in range(n_runs):
        ax.plot([x_a[i], x_b[i]],
                [data_a[i], data_b[i]],
                color='gray', linewidth=1, alpha=0.6,
                zorder=5, solid_capstyle='round')
    
    # C. Overlay means # TODO: is this necessary
    # mean_a = np.mean(data_a)
    # mean_b = np.mean(data_b)
    # ax.scatter([0, 1], [mean_a, mean_b], 
    #             color="#5DDDC6", edgecolor="#5DDDC6", s=90, marker='s', zorder=11, linewidth=1.8)

    # D. Significance annotation # TODO: improve so it does NOT cover the violins
    p_val = row['p_value_bonferroni']
    is_sig = p_val < 0.001
    
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    bracket_y = y_max - (y_range * 0.08)
    star_y = y_max - (y_range * 0.05)
    text_y = y_max - (y_range * 0.12)

    if is_sig:
        ax.plot([0, 1], [bracket_y, bracket_y], color='black', linewidth=1.4)
        ax.plot([0, 0], [bracket_y, bracket_y - (y_range*0.02)], color='black', linewidth=1.4)
        ax.plot([1, 1], [bracket_y, bracket_y - (y_range*0.02)], color='black', linewidth=1.4)
        
        ax.text(0.5, star_y, '***', ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')
        
        if p_val < 0.0001:
            p_str = f"{p_val:.1e}"
        else:
            p_str = f"{p_val:.4f}"
        
        ax.text(0.5, text_y, p_str, ha='center', va='top', fontsize=11, color='black', fontweight='bold')
    else:
        ax.text(0.5, bracket_y, 'ns', ha='center', va='bottom', fontsize=12, color='gray', style='italic')

    # E. Formatting
    ax.set_title(f"d = {row['cohen_d']:.2f}", fontsize=12, pad=12, fontweight='bold')
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([roi_a_name, roi_b_name], fontsize=12, weight='bold')
    ax.tick_params(axis='y', labelsize=9, direction='in', length=6)
    ax.tick_params(axis='x', length=0)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.spines['left'].set_linewidth(1.2)
    
    ax.grid(axis='y', linestyle='--', alpha=0.4, linewidth=0.8, zorder=0)

# Hide empty subplots
for j in range(n_tests, len(axes)):
    fig.delaxes(axes[j])

# *~*~*~*~ Optionally save # TODO: saving not GOOD
if save_fig:
    fig_name = f"sub-{subID:02d}_ses-{sessions}_block-{blocks}_space-{space}_contrast-{jobName}.png"
    fig_path = outPath / fig_name
    plt.savefig(fig_path, dpi = 300, bbox_inches = "tight")
    plt.close(fig)

# *~*~*~*~ Show
plt.show()