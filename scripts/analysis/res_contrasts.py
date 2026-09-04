#! /usr/bin/env python
# Time-stamp: <02-09-2026 m.utrosa@bcbl.eu>
"""
Conducts the 2nd level analysis on averaged contrast image

Returns:
- effect sizes, p-values, and t-tests for each region pair
- figure showing whether an ROI pair differs in its response to a specific contrast
"""

# Import python packages
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from itertools import combinations # generate all possible combos
                                   # order does not matter 

# Import custom-made functions
import config as c
from utils import extract_roi_array
from config import plotConf, apply_figure_style
apply_figure_style()

# MAIN TODOS: make it work for conditions regardless of the number of conditions
#             make it work for for averaging across voxels and runs
# TODO: 03b. could be extended to compare voxels (or voxel groups) per roi.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 00. Check that inputs are defined correctly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if c.average_voxels and c.average_runs:
    raise ValueError(
        "Statistical tests cannot be performed when averaging across BOTH runs and voxels. "
        "This results in a single scalar value per ROI (N=1). "
        "Please set either average_voxels=False or average_runs=False."
    )

if not c.average_voxels and not c.average_runs:
    raise ValueError(
        "To perform statistical tests we need one-dimensional arrays, "
        "which means that the selected values per contrast images "
        "have to be averaged EITHER across runs or voxels."
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Sum all contrast images into one (and save to disk)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
con_info = []
for sesID in c.sesIDs:
    for acqID in c.acqIDs:

        # Create a contrast dictionary
        con_dict = {}
        con_dict["sesID"] = sesID
        con_dict["acqID"] = acqID

        # Define the contrast file and add to dict
        con_name = c.con_filename
        con_fold = c.dataPath / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
        con_path = con_fold / con_name

        # Check that the file exists
        if con_path.exists():
            con_dict["con_path"] = con_path
        else:
            raise FileNotFoundError(f"The file does not exit: {con_path}")
        # Append the dict to a list of contrasts
        con_info.append(con_dict)

# Sum the contrasts images
con_img_sum = sum([nib.load(con["con_path"]).get_fdata() for con in con_info])
con_affine  = nib.load(con_info[0]["con_path"]).affine
print("\nAssuming all contrast images have the same affine.")

# Optionally, save the summed contrast
if c.save_summed:
    sum_name = f"sub-{c.subID:02d}_ses-{c.sessions}_acq-BLOCK{c.blocks}_summed-{c.con_filename}.gz"
    nib.save(
        nib.Nifti1Image(con_img_sum, con_affine),
        c.out_2nd / sum_name
        )
    print(f"\nSaved summed contrast image as {sum_name} to {c.out_2nd}.\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02a. Extract contrast values per voxel from each ROI.
# Shape of the extracted roi array: (n_runs, n_voxels)
# Size: {[[], []], # ROI 1: n_run lists with each sublist having n_voxel items
#        [[], []]}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initialize a dictionary to save extracted values
extracted_con = {name : [] for name in c.rois.keys()}
roi_names = list(extracted_con.keys())

# Iterate through each contrast image in the list of contrasts images
for con_dict in con_info:
        
        # Extract the subcortical arrays		
        mask_subcor, _, con__subcor_affine = extract_roi_array(
            c.subID,
            con_dict["sesID"],
            con_dict["acqID"],
            c.atlas_subcor_path,
            c.space,
            con_dict["con_path"],
            c.rois_subcortical,
            c.out_2nd,
            verbose=False,
            save=c.save_roi,
            average_voxels=False # Keeping this false for consistency
        )
        
        # Extract the cortical arrays
        mask_cor, _, con_cor_affine = extract_roi_array(
            c.subID,
            con_dict["sesID"],
            con_dict["acqID"],
            c.atlas_cor_path,
            c.space,
            con_dict["con_path"],
            c.rois_cortical,
            c.out_2nd,
            verbose=False,
            save=c.save_roi,
            average_voxels=False # Keeping this false for consistency
        )

        # Accumulate subcortical arrays
        for name in c.rois_subcortical.keys():
            extracted_con[name].append(mask_subcor[name])

        # Accumulate cortical arrays
        for name in c.rois_cortical.keys():
            extracted_con[name].append(mask_cor[name])

# Print shape of the raw extracted data
print("\n--- RAW EXTRACTED DATA ---")
for roi_name in roi_names:
    print(f"{roi_name}: {np.shape(extracted_con[roi_name])}")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02b. Transform the extracted contrast values.
# Average across runs or voxels: (n_runs, n_voxels)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
selected_cons = extracted_con

# Collapse runs and voxels:
if c.average_voxels and c.average_runs:
    run_voxel_average = {}
    for name, array_list in selected_cons.items():
        run_voxel_average[name] = np.mean(array_list) 
    selected_cons = run_voxel_average 

# Collapse runs: get a mean contrast value across runs
elif c.average_runs:
    run_average = {}
    for name, array_list in selected_cons.items():
        run_average[name] = np.mean(array_list, axis=0) 
    selected_cons = run_average

# Collapse voxels: get a mean contrast value across voxels
elif c.average_voxels:
    voxel_average = {}
    for name, array_list in selected_cons.items():
        voxel_average[name] = np.mean(array_list, axis=1) 
    selected_cons = voxel_average

# Print update on the structure of array
print(
    "\n--- TRANSFORMED EXTRACTED DATA ---",
    f"\nAveraged across runs: {c.average_runs}",
    f"\nAveraged across voxels: {c.average_voxels}\n")
for roi_name in roi_names:
    if np.isscalar(selected_cons[roi_name]):
        print(f"{roi_name}: {float(selected_cons[roi_name]):.4f}")
    else:
        print(f"{roi_name}: {np.shape(np.array(selected_cons[roi_name]))}")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03a. Compare extracted and transformed contrasts per ROI against ZERO
# RQ: Does this specific ROI distinguish between the defined contrast?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(
    "\nPerforming two-sided Wilcoxon rank tests per ROI against 0 "
    f"for contrast {c.con_name}.\n"
)
results_one_sample = []
for roi in roi_names:

    data = np.array(selected_cons[roi])

    # Ensure data is 1D
    if data.ndim == 0:
        raise ValueError(
            f"\nFor ROI '{roi}' the data is a scalar (N=1). "
            "Cannot perform statistics. Check averaging settings."
        )
    elif data.ndim > 1:
        raise ValueError(
            f"\nFor ROI '{roi}' the data is not 1D: {data.shape}. "
            "Check averaging parameters. Wilcoxon test requires one-dimensional input.")
    
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
        "roi"            : roi,
        "N"              : len(data),
        "mean_contrast"  : mean_val,
        "std_contrast"   : std_val,
        "stat"           : statistic,
        "p_value"        : p_value,
        "cohen_d"        : cohen_d
    })

# Create and display DataFrame
df_results_against0 = pd.DataFrame(results_one_sample)
print(df_results_against0.head())

# Calculate the number of tests for Bonferroni correction
n_tests = len(df_results_against0)
print(f"\nThe Bonferroni correction is applied for {n_tests} tests.")

# Adjust p-values with Bonferroni
df_results_against0['p_value_bonferroni'] = df_results_against0['p_value'] * n_tests

# Reduce the adjusted p-values that exceed 1 to 1
df_results_against0['p_value_bonferroni'] = df_results_against0['p_value_bonferroni'].clip(upper=1.0)

# Save dataframe
df_results_against0.to_csv(
    c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_con-{c.con_name}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-against0.csv",
    index=False
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03b. Compare pairs of ROIs to see if they respond differently
# RQ: Is the contrast estimate larger in one ROI than another?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(
    "\nPerforming two-sided, paired Wilcoxon rank tests "
    f"for contrast {c.con_name}.\n"
)

# Paired rank sum tests can only be done comparing arrays of equal length.
# Compare rois across runs because this ensures the same number of observations.
if c.average_voxels:
    results_paired = []
    for roi_a, roi_b in combinations(roi_names, 2):
        
        # Extract the paired arrays
        data_a = np.array(selected_cons[roi_a])
        data_b = np.array(selected_cons[roi_b])
        
        # Ensure lengths match
        if len(data_a) != len(data_b):
            continue

        # Ensure both data are 1D
        if data_a.ndim == 0 or data_b.ndim == 0:
            raise ValueError(
                f"\nFor ROI {roi_a} or {roi_b} the data is a scalar (N=1). "
                "Cannot perform statistics. Check averaging settings."
            )
        elif data_a.ndim > 1 or data_b.ndim > 1:
            raise ValueError(
                f"\nFor ROI {roi_a} or {roi_b} the data is not 1D: {data_a.shape} and {data_b.shape}."
                " Check averaging parameters. Wilcoxon test requires one-dimensional input.")
        
        # The null hypothesis: the two related paired samples (rois) come from the same distribution. 
        diff_data = np.round(data_a - data_b, decimals=8)
        res = wilcoxon(
            x = diff_data,             # Differences between two sets of measurements (1D)  
            y = None,                  # Keep undefined to avoid roundoff error 
            zero_method = "pratt",     # Conservative: includes zero-differences 
            correction  = False,       # Default
            alternative = "two-sided", # d = x 
            method      = "auto",      # Default
            axis        = 0,
            nan_policy  = "propagate"  # How to handle nan-values: propagate or omit 
        )

        # Calculate metrics
        mean_val  = np.mean(diff_data)
        std_val   = np.std(diff_data, ddof=c.ddof)
        cohen_val = mean_val / std_val if std_val > 0 else np.nan

        # Append metrics
        results_paired.append({
            "roi_a"   : roi_a,
            "roi_b"   : roi_b,
            "mean"    : mean_val,
            "std"     : std_val,
            "stat"    : res.statistic, # This is the sum of the ranks of the differences
                                       # above or below zero (for two-sided)
            "p_value" : res.pvalue,
            "cohen_d" : cohen_val
        })

    # Create and display the dataframe
    df_results_paired = pd.DataFrame(results_paired)
    print(df_results_paired.head())

    # Apply Bonferroni correction
    n_tests = len(df_results_paired)
    print(f"\nThe Bonferroni correction is applied for {n_tests} tests.")

    # Adjust p-vaues with Bonferroni correction
    df_results_paired['p_value_bonferroni'] = df_results_paired['p_value'] * n_tests

    # Reduce the adjusted p-values that exceed 1 to 1
    df_results_paired['p_value_bonferroni'] = df_results_paired['p_value_bonferroni'].clip(upper=1.0)

    # Save dataframe
    df_results_paired.to_csv(
        c.out_2nd / f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_con-{c.con_name}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_test-paired.csv",
        index=False
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04a. Plotting results from inferential tests in 03a
# RQ: Are the contrast values per ROI significantly different from zero?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print(
    "\nThe session and acquisition labels corresponding to each run are saved "
    f"in 'con_info' variable. The order of contrast images (filenames) in it "
    "corresponds to the order of runs in the data that is being plotted."
    )
color   = "#9EDDFB" 
df_plot = df_results_against0.copy()
n_rois  = len(df_plot)

# Grid
fig, axes = plt.subplots(
    int(np.ceil(n_rois / plotConf["cols"])),
    plotConf["cols"],
    figsize=plotConf["figsize"],
    constrained_layout=True,
    sharey=True
)
axes = axes.flatten()

# Get limits for y axis
all_values = []
for _, row in df_plot.iterrows():
    roi_name = row['roi']
    data = selected_cons[roi_name]
    all_values.extend(data)
global_min = min(all_values)
global_max = max(all_values)
range_val  = global_max - global_min
y_min_global = global_min - (range_val * 0.1)
y_max_global = global_max + (range_val * 0.1)

# Plotting loop
for idx, (_, row) in enumerate(df_plot.iterrows()):
    ax = axes[idx]
    roi_name = row['roi']
    
    # Get raw data
    data = selected_cons[roi_name]
    n_runs = len(data)
    
    # Create plotting dataframe
    plot_df = pd.DataFrame({
        'ROI': [roi_name] * n_runs,
        'Value': data,
        'Run': list(range(1, n_runs + 1))
    })

    # A. Plot violin
    sns.violinplot(
        data=plot_df, 
        x='ROI', 
        y='Value',
        ax=ax,
        color=color, 
        inner=None,
        linewidth=1.2,
        width=0.2, # Make the violin compact
        cut=0 # Limit the violin within the data range!
    )
    
    # B. Plot individual points: voxels or runs
    sns.swarmplot(
            data=plot_df,
            x='ROI',
            y='Value',
            hue='Run',
            ax=ax,
            size=6,
            zorder=10,
            palette="YlOrBr",
            legend=False
        )

    # C. Overlay mean, SEM bar, and zero-reference line
    mean_val = np.mean(data)
    sem_val = np.std(data, ddof=c.ddof) / np.sqrt(n_runs)
    
    # Mean is a large black dot with black edge
    ax.scatter(0.03, mean_val, 
               color="black", edgecolor='black', s=35, marker='o', 
               zorder=11, linewidth=0.6)
    
    # Plot error bar in orange
    ax.errorbar(0.03, mean_val, yerr=sem_val, 
                fmt='none', ecolor='#FF8C00', linewidth=1.2, 
                capsize=6, capthick=2, zorder=12)
    
    # Reference line at zero
    ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8, zorder=1)

    # D. Annotate significance (against zero) and add titles
    p_val  = row['p_value_bonferroni']
    is_sig = p_val < 0.05
    
    # Determine position for stars (above the max value or above the error bar)
    y_range = y_max_global - y_min_global
    
    # Place stars above the data
    star_y = y_max_global + (y_range * 0.05) 
    text_y = y_max_global + (y_range * 0.08)
    
    # Add text
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
        
        ax.text(0, star_y, stars, ha='center', va='center', 
                fontsize=plotConf["subplot_fontsize"], fontweight ="bold", color='black')
        ax.text(0, text_y, p_str, ha='center', va='bottom', 
                fontsize=plotConf["subplot_fontsize"], color='black')
    else:
        ax.text(0, star_y, "ns", ha='center', va='bottom', 
                fontsize=plotConf["subplot_fontsize"], color='gray', style='italic')
    
    # Title with mean contrast value and effect size
    ax.set_title(
        f"{roi_name} | μ={mean_val:.3f}, d={row['cohen_d']:.2f}",
        fontsize=plotConf["subplot_fontsize"],
        pad=10,
        va='top'
    )
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([])

    # Apply y limits
    ax.set_ylim(y_min_global, y_max_global)

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
    f'Contrast Estimate ({c.con_name})',
    va='center',
    rotation='vertical',
    fontsize=plotConf["fig_fontsize"]
)

# Optionally save and show the figure
if c.save_fig:
    fig_name = f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_con-{c.con_name}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}_one_sample.png"
    fig_path = c.out_2nd / fig_name
    plt.savefig(fig_path, dpi=plotConf["dpi"], bbox_inches="tight")
    plt.close(fig)

if c.show_fig:
    plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04b. Plotting
# RQ: Do mean contrast estimates differ for each ROI pair?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if c.average_voxels:
#     cb_palette = ["#FFBB28", "#8B84D8"]

#     # *~*~*~*~ Setup grid and axes
#     cols = plotConf["cols"]
#     rows = int(np.ceil(n_tests / cols))
#     fig, axes = plt.subplots(rows, cols, figsize=plotConf["figsize"], constrained_layout=True)
#     axes = axes.flatten()

#     # *~*~*~*~ Plotting loop
#     df_plot = df_results_paired.copy()
#     for idx, (_, row) in enumerate(df_plot.iterrows()):
#         ax = axes[idx]

#         # Get roi name
#         roi_a_name = row['roi_a']
#         roi_b_name = row['roi_b']
        
#         # Get raw data
#         data_a = selected_cons[roi_a_name]
#         data_b = selected_cons[roi_b_name]
#         n_runs = len(data_a)
        
#         # Create a dataframe for seaborn
#         # TODO: plot_df is probably wrong!!
#         plot_df = pd.DataFrame({
#             'ROI': [roi_a_name] * n_runs + [roi_b_name] * n_runs,
#             'Value': np.concatenate([data_a, data_b]),
#             'Run_ID': list(range(1, n_runs + 1)) * 2
#         })
        
#         # Ensure consistency in coloring
#         unique_rois = plot_df['ROI'].unique()
#         current_palette = {roi_a_name: cb_palette[0], roi_b_name: cb_palette[1]}

#         # A. Plot split violins
#         sns.violinplot(
#             data=plot_df, 
#             x='ROI', 
#             y='Value',
#             hue='ROI', 
#             ax=ax,
#             palette=current_palette, 
#             hue_order=[roi_a_name, roi_b_name],
#             split=True, 
#             inner=None,
#             linewidth=1.5,
#             alpha=0.8,
#             legend=False
#         )
        
#         # B. Plot individual points with controlled jitter + correct connections
#         rng = np.random.default_rng(42)
#         jitter_strength = 0.12

#         # deterministic jitter for pairing consistency
#         jitter_a = rng.uniform(-jitter_strength, jitter_strength, n_runs)
#         jitter_b = rng.uniform(-jitter_strength, jitter_strength, n_runs)

#         x_a = np.zeros(n_runs) + jitter_a
#         x_b = np.ones(n_runs) + jitter_b

#         # Scatter individual constrast estimate points
#         ax.scatter(x_a, data_a,
#                 color='black', s=35, alpha=0.8,
#                 edgecolor='black', linewidth=1, zorder=10)

#         ax.scatter(x_b, data_b,
#                 color='black', s=35, alpha=0.8,
#                 edgecolor='black', linewidth=1, zorder=10)

#         # Draw paired connections
#         for i in range(n_runs):
#             ax.plot([x_a[i], x_b[i]],
#                     [data_a[i], data_b[i]],
#                     color='gray', linewidth=1, alpha=0.6,
#                     zorder=5, solid_capstyle='round')
        
#         # C. Overlay means # TODO: is this necessary
#         # mean_a = np.mean(data_a)
#         # mean_b = np.mean(data_b)
#         # ax.scatter([0, 1], [mean_a, mean_b], 
#         #             color="#5DDDC6", edgecolor="#5DDDC6", s=90, marker='s', zorder=11, linewidth=1.8)

#         # D. Significance annotation # TODO: improve so it does NOT cover the violins
#         p_val = row['p_value_bonferroni']
#         is_sig = p_val < 0.001
        
#         y_min, y_max = ax.get_ylim()
#         y_range = y_max - y_min
#         bracket_y = y_max - (y_range * 0.08)
#         star_y = y_max - (y_range * 0.05)
#         text_y = y_max - (y_range * 0.12)

#         if is_sig:
#             ax.plot([0, 1], [bracket_y, bracket_y], color='black', linewidth=1.4)
#             ax.plot([0, 0], [bracket_y, bracket_y - (y_range*0.02)], color='black', linewidth=1.4)
#             ax.plot([1, 1], [bracket_y, bracket_y - (y_range*0.02)], color='black', linewidth=1.4)
            
#             ax.text(0.5, star_y, '***', ha='center', va='bottom', fontsize=plotConf["fig_fontsize"], fontweight='bold', color='black')
            
#             if p_val < 0.0001:
#                 p_str = f"{p_val:.1e}"
#             else:
#                 p_str = f"{p_val:.4f}"
            
#             ax.text(0.5, text_y, p_str, ha='center', va='top', fontsize=plotConf["subplot_fontsize"], color='black', fontweight='bold')
#         else:
#             ax.text(0.5, bracket_y, 'ns', ha='center', va='bottom', fontsize=plotConf["subplot_fontsize"], color='gray', style='italic')

#         # E. Formatting
#         ax.set_title(f"d = {row['cohen_d']:.2f}", fontsize=plotConf["subplot_fontsize"], pad=12, fontweight='bold')
#         ax.set_xlabel("")
#         ax.set_ylabel("")
#         ax.set_xticks([0, 1])
#         ax.set_xticklabels([roi_a_name, roi_b_name], fontsize=plotConf["subplot_fontsize"], weight='bold')
#         ax.tick_params(axis='y', labelsize=9, direction='in', length=6)
#         ax.tick_params(axis='x', length=0)
        
#         ax.spines['top'].set_visible(False)
#         ax.spines['right'].set_visible(False)
#         ax.spines['bottom'].set_linewidth(1.2)
#         ax.spines['left'].set_linewidth(1.2)
        
#         ax.grid(axis='y', linestyle='--', alpha=0.4, linewidth=0.8, zorder=0)

#     # Hide empty subplots
#     for j in range(n_tests, len(axes)):
#         fig.delaxes(axes[j])

#     # *~*~*~*~ Optionally save # TODO: saving not GOOD
#     if c.save_fig:
#         fig_name = f"sub-{c.subID:02d}_ses-{c.sessions}_block-{c.blocks}_space-{c.space}_job-{c.jobName}_con-{c.con_name}_avgVox-{c.average_voxels}_avgRun-{c.average_runs}.png"
#         fig_path = c.out_2nd / fig_name
#         plt.savefig(fig_path, dpi=plotConf["dpi"], bbox_inches = "tight")
#         plt.close(fig)

#     if c.show_fig:
#         plt.show()