#!/usr/bin/env python
# Time-stamp: <15-05-2026 m.utrosa@bcbl.eu>
# conda activate localizer_fMRI
"""
Load participant's logs to inspect difficulty of frequency counting task.
"""

# ---------- PREP
# Python packages
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

# Customized modules
from scripts import grabber
import percentCorr as pc

# ---------- PATHS
# Get logfiles
homePath = Path("/home/mutrosa/Documents/projects")
logPath = homePath / "devLoc" / "data_logs" 
designPath = homePath / "auditory_paradigms" / "detection_accuracy" / "selected_trials"
plotsDir = homePath / "devLoc" / "data_logs" / "plots"

# ---------- SESSION INFO
no_blocks = int(input("Enter the number of blocks:"))
comboID   = int(input("Enter the index of exp_parameter_combo.csv:"))
sesID = int(input("Enter the session ID:"))
subID = int(input("Enter the subject ID:"))

# ---------- EXTRACT DATA
df_freqDiff = pd.DataFrame()
df_devNum = pd.DataFrame()

for i in range(no_blocks):
	block_idx = i + 1

	# Get responses and frequency difference by design
	freq_dev_no_design = pc.extract_design_data(
	    trials_dir = designPath,
	    combo = comboID,
	    block = block_idx,
	    target_column = "freq_dev_no")

	freq_diff_design   = pc.extract_design_data(
	    trials_dir = designPath,
	    combo = comboID,
	    block = block_idx,
	    target_column = "freq_diff")

	# Get responses by participant
	freq_dev_no_participant = pc.extract_participant_data(
	    project_dir = logPath / f"bids_output_{comboID:003d}",
	    subject = subID,
	    session = sesID,
	    task = "freqDev",
	    block = block_idx)

	# 01. Calculate percentage correct per block
	# Make sure that the data is of the same type
	resp_design_int  = freq_dev_no_design.astype(int)
	resp_design_str  = resp_design_int.astype(str)
	resp_subject_str = freq_dev_no_participant.astype(str)

	# Count no. of correct responses
	resp_bool = resp_design_str == resp_subject_str
	resp_type = []
	for r in resp_bool:
		if r == False:
			resp_type.append("Wrong")
		else:
			resp_type.append("Correct")

	# 02. Extract frequency differences counts per block
	flat_freq_diff = freq_diff_design.explode()
	freqDiff_counts = flat_freq_diff.groupby(flat_freq_diff).count()
	freqDiff_dict = freqDiff_counts.to_dict()

	# 03. Extract frequency deviants number counts per block
	flat_freq_dev_no = freq_dev_no_design.explode()
	freqDevNo_counts = flat_freq_dev_no.groupby(flat_freq_dev_no).count()
	devNum_dict = freqDevNo_counts.to_dict()

	# Combine
	devNum = pd.concat([pd.Series(resp_type, dtype="str", name="response"), freq_dev_no_design], axis=1)
	devNum["block_idx"] = block_idx

	freqDiff = pd.concat([pd.Series(resp_type, dtype="str", name="response"), freq_diff_design], axis=1)
	freqDiff["block_idx"] = block_idx

	df_freqDiff = pd.concat([df_freqDiff, freqDiff])
	df_devNum = pd.concat([df_devNum, devNum])

# ---------- RESHAPE THE DATA FOR PLOTTING
# Explode data for frequency difference
df_long_freqDiff = df_freqDiff.explode("freq_diff")

# Ensure correct data types
df_devNum["freq_dev_no"] = df_devNum["freq_dev_no"].astype(int)
df_devNum["block_idx"]   = df_devNum['block_idx'].astype(int)
df_devNum["response"] = df_devNum["response"].astype(str)

df_long_freqDiff['freq_diff'] = df_long_freqDiff['freq_diff'].astype(int)
df_long_freqDiff['block_idx'] = df_long_freqDiff['block_idx'].astype(int)
df_long_freqDiff['response']  = df_long_freqDiff['response'].astype(str)

# Get unique blocks and sort them
blocks = sorted(df_long_freqDiff['block_idx'].unique())
n_blocks = len(blocks)

# Determine grid layout (e.g., 3 columns, dynamic rows)
n_cols = 2
n_rows = int(np.ceil(n_blocks / n_cols))

# ---------- PLOTTING 01: COUNT OF TRIALS WITH CORRECT/WRONG RESPONSES PER FREQ DIFF
plot_data = df_long_freqDiff.groupby(['block_idx', 'freq_diff', 'response']).size().reset_index(name='count')

# Pallete and theme
palette = {'Correct': '#00CED1', 'Wrong': '#FFD700'} # Cyan and Yellow
sns.set_theme(style="white")

# Figure with barplot subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
axes = axes.flatten()
for i, block in enumerate(blocks):
    ax = axes[i]
    
    # Filter data for current block
    block_data = plot_data[plot_data['block_idx'] == block]

    sns.barplot(
        data=block_data,
        x='freq_diff',
        y='count',
        hue='response',
        hue_order=["Correct", "Wrong"],
        palette=palette,
        ax=ax
    )

    # Styling the subplot
    ax.set_title(f'Block {block}', fontsize=14, weight='bold', pad=10)
    ax.set_xlabel('Difference Between Deviant and Standard Frequency', fontsize=12, labelpad=5)
    ax.set_ylabel('Count', fontsize=12, labelpad=5)
    
    # Add grid lines for easier reading
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    ax.set_axisbelow(True)

    legend = ax.get_legend()
    if legend:
        legend.remove()

# Create a single global legend
handles, labels = ax.get_legend_handles_labels() # Get from last active ax
# If the last ax was empty, we might need to create dummy handles
if not handles:
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=palette['Correct'], edgecolor='black'), 
               Patch(facecolor=palette['Wrong'], edgecolor='black')]
    labels = ['Correct', 'Wrong']

fig.legend(handles, labels, title='Response', 
           loc='upper right', bbox_to_anchor=(0.4, 0.46, 0.5, 0.5), 
           ncol=2, frameon=True, fontsize=12, title_fontsize=14)

# Global Title
fig.suptitle(
    f'Response Accuracy by Frequency Difference per Block\n(sub-{subID:02d}, ses-{sesID:02d}, trials-{comboID:003d})', 
    fontsize=18, weight='bold', y=0.98
)

plt.tight_layout(rect=[0, 0, 1, 0.96])
outPath = plotsDir / f"sub-{subID:02d}_ses-{sesID:02d}_trials-{comboID:003d}_freqDiff_bar.png"
plt.savefig(outPath, dpi=300, bbox_inches="tight")
plt.show()

# ---------- PLOTTING 02: COUNT OF TRIALS WITH CORRECT/WRONG RESPONSES PER FREQ DEV NO
plot_data = df_devNum.groupby(['block_idx', 'freq_dev_no', 'response']).size().reset_index(name='count')

# Pallete and theme
palette = {'Correct': '#00CED1', 'Wrong': '#FFD700'} # Cyan and Yellow
sns.set_theme(style="white")

# Figure with barplot subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
axes = axes.flatten()
for i, block in enumerate(blocks):
    ax = axes[i]
    
    # Filter data for current block
    block_data = plot_data[plot_data['block_idx'] == block]

    sns.barplot(
        data=block_data,
        x='freq_dev_no',
        y='count',
        hue='response',
        hue_order=["Correct", "Wrong"],
        palette=palette,
        ax=ax
    )

    # Styling the subplot
    ax.set_title(f'Block {block}', fontsize=14, weight='bold', pad=10)
    ax.set_xlabel('Number of Frequency Deviants per Trial', fontsize=12, labelpad=5)
    ax.set_ylabel('Count', fontsize=12, labelpad=5)
    
    # Add grid lines for easier reading
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    ax.set_axisbelow(True)

    legend = ax.get_legend()
    if legend:
        legend.remove()

# Create a single global legend
handles, labels = ax.get_legend_handles_labels()

# If the last ax was empty, we might need to create dummy handles
if not handles:
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=palette['Correct'], edgecolor='black'), 
               Patch(facecolor=palette['Wrong'], edgecolor='black')]
    labels = ['Correct', 'Wrong']

fig.legend(handles, labels, title='Response', 
           loc='upper right', bbox_to_anchor=(0.4, 0.46, 0.5, 0.5), 
           ncol=2, frameon=True, fontsize=12, title_fontsize=14)

# Global Title
fig.suptitle(
    f'Response Accuracy by Number of Frequency Deviants per Block\n(sub-{subID:02d}, ses-{sesID:02d}, trials-{comboID:003d})', 
    fontsize=18, weight='bold', y=0.98
)

plt.tight_layout(rect=[0, 0, 1, 0.96])
outPath = plotsDir / f"sub-{subID:02d}_ses-{sesID:02d}_trials-{comboID:003d}_devNum_bar.png"
plt.savefig(outPath, dpi=300, bbox_inches="tight")
plt.show()