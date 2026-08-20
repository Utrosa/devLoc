# Inspect motion outliers
import json
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

block = 2
sesID = 2
subID = 5
task = "timDev"
acqID = f"BLOCK{block}"
denosing = True
save = True

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

# 00. Files and paths
# Path to the project directory
homePath = Path('/home/mutrosa/mutrosa/Documents/projects/devLoc')

# Get info about the filtered rigid body movements: trasn & rot
mriFold = homePath / "data_MRI" / 'derivatives' / f"NORDIC-{denosing}" / "derivatives"
mriFunc = mriFold / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / "func"
mriName = f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_desc-confounds_timeseries.json"
mriPath = mriFunc / mriName

# Get the filtered rigid body movements: trasn & rot
artFold = homePath / "data_physio" / "artifacts" / f"NORDIC-{denosing}"
artName = f"sub-{subID:02d}_ses-{sesID:02d}_acq-BLOCK{block}_movpar.txt"
artPath = artFold / artName

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

# Translations
motion[["trans_x", "trans_y", "trans_z"]].plot(ax=axes[0])
axes[0].set_ylabel("Translation [mm]")
axes[0].set_title("Motion Parameters")

# Rotations
motion[["rot_x", "rot_y", "rot_z"]].plot(ax=axes[1])
axes[1].set_ylabel("Rotation [rad]")
axes[1].set_xlabel("Volume")

plt.tight_layout()
plt.show()