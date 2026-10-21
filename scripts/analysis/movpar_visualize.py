#! /usr/bin/env python
# Time-stamp: <20-10-2026 m.utrosa@bcbl.eu>
"""
Visualizes extent of motion per block and session.
Motion is estimated with rigid-body transformations from fMRIPrep.
These plots could be useful for:
* determining individual rot & trans thresholds for outlier detection.
* comparing with plots from Nipype's ArtifactDetection() algorithm.

Possible adaptations:
* y axis limits
* csf & wm confounds

https://github.com/nipy/nipype/blob/master/nipype/algorithms/rapidart.py
"""
import json
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

blocks = [1, 2, 3, 4]
sessions = [2, 3, 4, 5, 6, 7]
subID = 5
task = "timDev"
denosing = True
save = True
show = False

# WARNING!
# Ensure that the order of the confounds keys is the same as in the list specified in the
# filter artifacts function! Because the ".movpar" file has no headers, we do not know which
# column is which and pandas will filter the original "timeseries.tsv" given the list of keys.
keys = [
    'rot_x',
    'rot_y',
    'rot_z'
    'trans_x', 
	'trans_y',
	'trans_z'
    ]
# csf_wm_keys = [
# 			"csf",
# 			"csf_derivative1",
# 			"csf_derivative1_power2", 
# 			"csf_power2", "white_matter",
# 			"white_matter_derivative1",
# 			"white_matter_derivative1_power2",
# 			"white_matter_power2",
# 			"csf_wm"
# 		]

# ------------------------------------------------------------------------------------------------------------------
# BELOW: DO NOT MODIFY
# ------------------------------------------------------------------------------------------------------------------
for sesID in sessions:
    for block in blocks:
        acqID = f"BLOCK{block}"
        
        # 00. Files and paths
        # Path to the project directory
        # homePath = Path('/home/mutrosa/mutrosa/Documents/projects/devLoc')
        homePath = Path('/home/mutrosa/Documents/projects/devLoc')

        # Get info about the filtered rigid body movements: trasn & rot
        mriFold = homePath / "data_MRI" / 'derivatives' / f"NORDIC-{denosing}" / "derivatives"
        mriFunc = mriFold / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / "func"
        mriName = f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_desc-confounds_timeseries.json"
        mriPath = mriFunc / mriName

        # Get the filtered rigid body movements: trasn & rot
        artFold = homePath / "data_physio" / "artifacts" / f"NORDIC-{denosing}"
        artName = f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-BLOCK{block}_movpar.txt"
        artPath = artFold / artName
        outPath = artFold / "movpar_plots"
        outPath.mkdir(parents=True, exist_ok=True)

        # 01. Load and read
        # Load motion parameter file
        motion = pd.read_csv(
            artPath,
            sep=r"\s+",
            header=None,
            names=["trans-x", "trans-y", "trans-z", "rot-x", "rot-y", "rot-z"]
        )

        # 02. Check order and unit of parameters
        with open(mriPath, "r") as f:
            confound_meta = json.load(f)

        for key in confound_meta:
            if "trans" in key or "rot" in key:
                print(f"\n{key}")
                print(confound_meta[key])

        # 03. Plot the movement
        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # Translations (mm)
        motion[["trans-x", "trans-y", "trans-z"]].plot(ax=axes[0])
        axes[0].set_ylabel("Translation [mm]")
        axes[0].set_title(f"sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-BLOCK{block}")
        axes[0].set_ylim(-0.35, 0.35)
        
        # Add lines at an initial default threshold of 0.3
        axes[0].axhline(0.1, color='black', linestyle='--', alpha=0.8, zorder=1)
        axes[0].axhline(-0.1, color='black', linestyle='--', alpha=0.8, zorder=1)

        axes[0].axhline(0.04, color='red', linestyle='--', alpha=0.8, zorder=1)
        axes[0].axhline(-0.04, color='red', linestyle='--', alpha=0.8, zorder=1)

        # Rotations (degrees)
        motion[["rot-x", "rot-y", "rot-z"]].plot(ax=axes[1])
        axes[1].set_ylabel("Rotation [rad]")
        axes[1].set_xlabel("Volume")
        axes[1].set_ylim(-3, 3)

        # Add lines at an initial default threshold of 0.3
        axes[1].axhline(0.3, color='black', linestyle='--', alpha=0.8, zorder=1)
        axes[1].axhline(-0.3, color='black', linestyle='--', alpha=0.8, zorder=1)

        axes[1].axhline(0.7, color='red', linestyle='--', alpha=0.8, zorder=1)
        axes[1].axhline(-0.7, color='red', linestyle='--', alpha=0.8, zorder=1)

        # Style
        plt.tight_layout()

        # Save and show the plots
        if save:
            plt.savefig(outPath / f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-BLOCK{block}_movpar_plot.png", dpi=300, bbox_inches="tight")
        if show:
            plt.show()