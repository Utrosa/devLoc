#!/usr/bin/env python
# Time-stamp: <15-05-2026 m.utrosa@bcbl.eu>
# conda activate localizer_fMRI
"""
Load participant's logs to inspect difficulty of frequency counting task.
Wrong plotting.
"""

# ---------- PREP
# Python packages
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
freqDiff_data = []
devNum_data = []
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
	no_same = resp_design_str == resp_subject_str
	no_corr = no_same.sum() # True values are summed
	no_all  = len(freq_dev_no_design)
	perCorr_block = no_corr / no_all

	# 02. Extract frequency differences counts per block
	flat_freq_diff = freq_diff_design.explode()
	freqDiff_counts = flat_freq_diff.groupby(flat_freq_diff).count()
	freqDiff_dict = freqDiff_counts.to_dict()

	# 03. Extract frequency deviants number counts per block
	flat_freq_dev_no = freq_dev_no_design.explode()
	freqDevNo_counts = flat_freq_dev_no.groupby(flat_freq_dev_no).count()
	devNum_dict = freqDevNo_counts.to_dict()

	# Save the extracted data
	row_data = {
		"block_idx" : block_idx,
		"perCorr" : perCorr_block
	}
	row_data.update(freqDiff_dict)
	freqDiff_data.append(row_data)

	row_data = {
		"block_idx" : block_idx,
		"perCorr" : perCorr_block
	}
	row_data.update(devNum_dict)
	devNum_data.append(row_data)

df_freqDiff = pd.DataFrame(freqDiff_data)
df_devNum = pd.DataFrame(devNum_data)

# ---------- RESHAPE THE DATA FOR PLOTTING
# Identify the columns that are frequency differences (everything except block_idx and perCorr)
meta_cols = ['block_idx', 'perCorr']
freq_diff_cols = [col for col in df_freqDiff.columns if col not in meta_cols]
dev_num_cols = [col for col in df_devNum.columns if col not in meta_cols]

# Melt the dataframe
df_long_freqDiff = df_block.melt(
    id_vars=meta_cols, 
    value_vars=freq_diff_cols, 
    var_name='freq_diff', 
    value_name='count'
)
df_long_devNum   = df_block.melt(
    id_vars=meta_cols, 
    value_vars=dev_num_cols, 
    var_name='freq_dev_no', 
    value_name='count'
)

# ---------- PLOTTING 01
sns.set_theme(style="white")
plt.figure(figsize=(12, 7))

scatter = sns.scatterplot(
    data=df_long_freqDiff, 
    x='freq_diff', 
    y='perCorr', 
    size='count',
    sizes=(50, 500), # Min and max bubble size (tuple)
    hue='freq_diff',
    palette='Set2',
    alpha=0.6, 
    edgecolor='black',
    linewidth=0.5,
    legend=False
)

# Labels and text
plt.axhline(0.5, color='red', linestyle='--', alpha=0.5, label='Chance Level')
plt.title(f'sub-{subID:02d}, ses-{sesID:02d}, trials-{comboID:003d}', fontsize=20, pad=20, weight='bold')
plt.xlabel('Difference Between Standard and Deviant Frequency', fontsize=18, labelpad=15)
plt.ylabel('Percentage of Correct Responses', fontsize=18, labelpad=15)

# Ticks
unique_diffs = sorted(df_long_freqDiff['freq_diff'].unique())
plt.xticks(unique_diffs)

plt.tight_layout()

outPath = plotsDir / f"sub-{subID:02d}_ses-{sesID:02d}_trials-{comboID:003d}_freqDiff.jpg"
plt.savefig(
	outPath,
	dpi=300,
	bbox_inches="tight"
	)
plt.show()

# ---------- PLOTTING: 02
sns.set_theme(style="white")
plt.figure(figsize=(12, 7))

scatter = sns.scatterplot(
    data=df_long_devNum, 
    x='freq_dev_no', 
    y='perCorr', 
    size='count',
    sizes=(50, 500), # Min and max bubble size (tuple)
    hue='freq_dev_no',
    palette='Set2',
    alpha=0.6, 
    edgecolor='black',
    linewidth=0.5,
    legend=False
)

# Labels and text
plt.axhline(0.5, color='red', linestyle='--', alpha=0.5, label='Chance Level')
plt.title(f'sub-{subID:02d}, ses-{sesID:02d}, trials-{comboID:003d}', fontsize=20, pad=20, weight='bold')
plt.xlabel('The number of frequency deviants per trial', fontsize=18, labelpad=15)
plt.ylabel('Percentage of Correct Responses', fontsize=18, labelpad=15)

# Ticks
unique_diffs = sorted(df_long_devNum['freq_dev_no'].unique())
plt.xticks(unique_diffs)

plt.tight_layout()

outPath = plotsDir / f"sub-{subID:02d}_ses-{sesID:02d}_trials-{comboID:003d}_devNum.jpg"
plt.savefig(
	outPath,
	dpi=300,
	bbox_inches="tight"
	)
plt.show()