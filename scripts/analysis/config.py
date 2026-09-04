#! /usr/bin/env python
# Time-stamp: <04-09-2026 m.utrosa@bcbl.eu>
"""
Configuration for the following scripts:
- resample_atlas.py
- res_contrasts.py
"""
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Activate python environment and import packages
# Citrix: source activate nipypee
# Local : conda activate nipypee
import yaml
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Specify type of denoising in preproc, the 1st level analysis, 
#     and the conditions modelled
# J1: whenwhat (timDev vs freqDev) => no specialization to scale
# => suboptimal as it only captures responses to wide patterns (grouping)
# J2: when11where (abs timDev vs freqDev)
# => captures specialization to temporal scale
denoising = True # NORDIC True or False
jobName   = "whenwhat" # when11where

# Conditions have to be in the order of beta images
# Please check the names in the SPM design matrix
if jobName == "whenwhat":
    conditions = ["timDev", "freqDev"]
    conditions_int = [1, 2]
elif jobName == "when11where":
    conditions = [4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125]
    conditions_int = list(range(1,12))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03. Specify data handling and plotting preferences
save_roi       = False  # applies to extracted ROI arrays
show_fig       = True   # applies to figures with statistical results
save_fig       = True   # applies to figures with statistical results
save_summed    = False  # applies to the contrasts images
save_averaged  = False  # averaged beta arrays
average_voxels = True  # CONTRASTS: If True, one value (array) across VOXELS.
                        # If both are False, the extracted roi array has shape (n_runs, n_voxels)
average_runs   = False   # If True, collapse runs and return a mean across runs.                  
remove_empty   = False  # Remove or not empty arrays (e.g.: If we do not average across voxels, 
					    # do we, when averaging across runs, include voxels that have zero 
					    # beta values or not?)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04. Specify experiment info
subID  = 5
anatID = 2
space  = "T1w" #TODO: What is the differences between T1w and T1wFOV
task   = "timDev"
sesIDs = [2, 3, 4, 5, 6, 7] # 2, 3, 4, 5, 6, 7
sessions = 234567 # appears in the filenames
acqIDs = ["BLOCK1", "BLOCK2", "BLOCK3", "BLOCK4"] # "FUNLOC" 
blocks = "1234" # appears in the filenames

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 05. Specify project directories
# homePath  = Path("/home/mutrosa/mutrosa/Documents/projects/devLoc") # Citrix
homePath  = Path("/home/mutrosa/Documents/projects/devLoc")           # Local
dataDir   = homePath / "results" # path to 1st level analysis results
dataPath  = dataDir / jobName / f"NORDIC-{denoising}" / "1stLevel"
outDir    = homePath / "tests" / jobName / f"NORDIC-{denoising}"
out_2nd   = outDir / "2ndLevel"
out_1st   = outDir / "1stLevel" / "visualization"
spmt_out  = homePath / "results" / "visualization"

# Create missing output directories
outDir.mkdir(parents=True, exist_ok=True)
out_2nd.mkdir(parents=True, exist_ok=True)
out_1st.mkdir(parents=True, exist_ok=True)
spmt_out.mkdir(parents=True, exist_ok=True)

# Raw data paths
mriPath  = homePath / "data_MRI" / "derivatives" / f"NORDIC-{denoising}" / "derivatives" 
anatPath = mriPath / f"sub-{subID:02d}" / f"ses-{anatID:02d}" / "anat"
funcPath = mriPath / f"sub-{subID:02d}"
freesurfer_dir = mriPath / "sourcedata" / "freesurfer" / f"sub-{subID:02d}_ses-{anatID:02d}" / "mri"

# Filenames and folders of the 1st level analysis output
# The 1st level results have to be resampled prior to visualization
con_name = "timDev-freqDev"                 # contrast label
con_filename  = "con_space-T1wFOV_0001.nii" # image
beta_filename = "beta_space-T1wFOV"         # image
spm_filename  = ""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 06. Specify atlas and roi info
with open("rois.yaml", "r") as f:
    data = yaml.safe_load(f)
rois_cortical    = data["cortical"]
rois_subcortical = data["subcortical"]
rois = data["cortical"] | data["subcortical"]

# Get Sitek's subcortical atlas
atlas_subcor_name = f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
atlas_subcor_path = homePath / "templates" / atlas_subcor_name

# Get FreeSurfer's parcellation: Destrieux Atlas
atlas_cor_name  = f"aparc.a2009s+aseg_NORDIC-{denoising}_space-{space}.nii.gz"
atlasPath       = homePath / "templates"
atlas_cor_path  = atlasPath / atlas_cor_name

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 07. Specify plotting style settings
def apply_figure_style():
    sns.set_context("paper", font_scale=1.3)
    sns.set_style("white")
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
    plt.rcParams['axes.linewidth'] = 1.2
plotConf = {
    "cols"             : 2,
    "figsize"          : (12, 15),
    "dpi"              : 300,
    "fig_fontsize"     : 14,
    "subplot_fontsize" : 12
}
plot_rois   = ["IC-L", "IC-R", "MGB-L", "MGB-R", "A1-L", "A1-R"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 08. Specify statistical tests
ddof = 1 # 1 = sample SD (with Bessel’s correction); 0 = population SD