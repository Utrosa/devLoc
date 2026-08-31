#! /usr/bin/env python
# Time-stamp: <31-08-2026 m.utrosa@bcbl.eu>
"""
Configuration for the following scripts:
- resample_atlas
- contrasts
- spmT
- betas
"""
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Citrix: source activate nipypee
# Local : conda activate nipypee

# Import python packages
import bids
import yaml
import warnings
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from nilearn.image import resample_to_img # for atlas img
from nibabel.orientations import aff2axcodes, io_orientation # for atlas img

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

# Conditions have to be in the order of beta images
# Please see names in the SPM design matrix
if jobName == "whenwhat":
    conditions = ["timDev", "freqDev"]
    conditions_int = [1, 2]
elif jobName == "when11where":
    conditions = [4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125]
    conditions_int = list(range(1,12))

# *~*~*~*~ Summation and plotting preferences *~*~*~*~
save_roi       = False # applies to extracted ROI arrays
save_fig       = True  # applies to figures with statistical results
save_average   = False # Applies to the summed contrast or beta arrays
average_voxels = False # Do we average voxels per ROI or not? If True, the code is wrong.
                       # TODO: Fix summation over ROIs. # TODO: Make it work for False for betas

# Only relevant for betas
remove_empty   = False # Remove or not empty arrays (e.g.: If we do not average across voxels, 
					   # do we, when averaging across runs, include voxels that have zero 
					   # beta values or not?)

# *~*~*~*~ Project directories *~*~*~*~
homePath  = Path("/home/mutrosa/mutrosa/Documents/projects/devLoc")
dataDir   = homePath / "results"
outDir    = homePath / "tests"
dataPath  = dataDir / jobName / f"NORDIC-{denoising}" / "1stLevel"
data_1stLevel = dataDir / jobName / f"NORDIC-{denoising}" / "1stLevel"

# !FOR CONTRASTS
 = outDir / jobName / f"NORDIC-{denoising}" / "2ndLevel"
outPath_con.mkdir(parents=True, exist_ok=True)
contrast_label = "con_space-T1wFOV_0001.nii" # name of the files from 1st Level Analysis

# !FOR BETAS
out_dir_beta  = outDir / jobName / f"NORDIC-{denoising}" / "1stLevel" / "visualization"
out_dir_beta.mkdir(parents=True, exist_ok=True)

# *~*~*~*~ Experiment info *~*~*~*~
subID  = 5
anatID = 2
space  = "T1w" #TODO: What is the differences between T1w and T1wFOV
task   = "timDev"
sesIDs = [2, 3, 4, 5, 6, 7] # 2, 3, 4, 5, 6, 7
sessions = 234567 # appears in the filenames
acqIDs = ["BLOCK1", "BLOCK2", "BLOCK3", "BLOCK4"] # "FUNLOC" 
blocks = "1234" # appears in the filenames
plot_rois = ["IC-L", "IC-R", "MGB-L", "MGB-R", "A1-L", "A1-R"]

# Raw data paths
mriPath  = homePath / "data_MRI" / "derivatives" / f"NORDIC-{denoising}" / "derivatives" 
anatPath = mriPath / f"sub-{subID:02d}" / f"ses-{anatID:02d}" / "anat"
funcPath = mriPath / f"sub-{subID:02d}"
freesurfer_dir = mriPath / "sourcedata" / "freesurfer" / f"sub-{subID:02d}_ses-{anatID:02d}" / "mri"

# *~*~*~*~ Atlases *~*~*~*~
with open("rois.yaml", "r") as f:
    data = yaml.safe_load(f)
rois_cortical    = data["cortical"]
rois_subcortical = data["subcortical"]
rois = data["subcortical"] | data["cortical"]

# Get Sitek's subcortical atlas
atlas_subcor_name = f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
atlas_subcor_path = homePath / "templates" / atlas_subcor_name

# Get FreeSurfer's parcellation: Destrieux Atlas
atlas_cor_name  = f"aparc.a2009s+aseg_NORDIC-{denoising}_space-{space}.nii.gz"
atlasPath       = homePath / "templates"
atlas_cor_path  = atlasPath / atlas_cor_name

# !FOR SPM t-images
atlas_path   = homePath / "templates" / f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
out_dir_spmt = homePath / "results" / "visualization"
out_dir_spmt.mkdir(exist_ok=True)