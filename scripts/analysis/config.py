#! /usr/bin/env python
# Time-stamp: <10-09-2026 m.utrosa@bcbl.eu>
"""
Configuration for the following scripts:
- resample_atlas.py
- resample_outputs.py
- res_contrasts.py
- res_betas.py
- res_spmT.py
- analysis_task-timDev.py
"""
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Activate python environment and import packages
# Citrix: source activate nipypee
# Local : conda activate nipype
import yaml
import seaborn as sns
from pathlib import Path
from utils import get_base_dirs
import matplotlib.pyplot as plt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Define the pipeline 
develop_mode = True  # developping mode
denoising    = True  # denoising from prepro (NORDIC)
verbose      = False

# How do we model deviant events? See deviants.yaml
timDev_jobName  = "when22"
freqDev_jobName = "what" # False or "what"
if freqDev_jobName:
    jobName  = timDev_jobName + freqDev_jobName
else:
    jobName = timDev_jobName

# Conditions determine the beta order. Check consistency with SPM design matrix.
with open("deviants.yaml", "r") as f:
    devs = yaml.safe_load(f)
timDev   = devs["timDev"][timDev_jobName]
if freqDev_jobName:
    freqDev  = devs["freqDev"][freqDev_jobName]
    conditions = timDev["conditions"] + freqDev["conditions"]
    contrast_weights = timDev["weights"] + freqDev["weights"]
else:
    conditions = timDev["conditions"]
    contrast_weights = timDev["weights"]

# Check correctness
conditions_int = list(range(1, len(conditions) + 1))
if not len(conditions_int) == len(conditions):
    raise ValueError(
        "The conditions do not match in length. Check configuration/yaml.")

# Design parameters for timing deviancy regressors
absolute = timDev["absolute"] # Absolute or nomibal timing deviancy regressors?
binary   = timDev["binary"]   # If True, magnitude (size) of timing deviants is not considered.
groups   = timDev["groups"]   # If False, no grouping. Values must be integers.

# Check if values are of correct type
if groups:
    for group_name, group_values in groups.items():
        for i, val in enumerate(group_values):
            if not isinstance(val, int):
                raise TypeError(
                    f"Error in group '{group_name}': Value at index {i} is {val!r} "
                    f"(type: {type(val).__name__}), expected int. "
                    "All values must be integers for bisect to work correctly."
                )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03a. Specify physiological regressors options
biopac = 1 # excludes (0) or includes (1) BIOPAC physiological regressorsx

# Confounds to be filtered from fMRIPrep
# None defaults to rot, trans, csf, and wm.
# The rigid body keys must be in order in which FSL expects them
# https://fsl.fmrib.ox.ac.uk/fsl/docs/registration/mcflirt.html
# confound_keys = [ 
#     "csf", 
#     "csf_derivative1", 
#     "csf_derivative1_power2", 
#     "csf_power2", 
#     "white_matter", 
#     "white_matter_derivative1", 
#     "white_matter_derivative1_power2", 
#     "white_matter_power2", 
#     "csf_wm"
# ]
confound_keys = ['rot_x', 'rot_y', 'rot_z', 'trans_x', 'trans_y', 'trans_z']

# Physiological regressors
tapas_cols = [f"RETROICOR_Cardiac_{i+1}" for i in range(6)] + \
             [f"RETROICOR_Respiratory_{i+1}" for i in range(8)] + \
             [f"RETROICOR_Multiplicative_{i+1}" for i in range(4)]

# Rapidart nipype node for motion artifact detection
artDetect = True
if artDetect:
    zintensity_thresh = 3   # detect images that deviate from the mean
    rot_thresh        = 0.3 # in radians
    trans_thresh      = 0.3 # in mm

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03b. Specify 1st level analysis options
concat    = True  # If False, treats runs as a single continuous series
hrf_dervs = [0, 0]
volterra  = False
smoothing = None  # Set the Gaussian filter width in mm 2.5; defaults to None

# Contrast specification
contrast  = True
contrasts = [(jobName, 'T', conditions, contrast_weights)]

# 03c. Specify data handling and plotting preferences for 1st level results
save_roi       = False  # applies to extracted ROI arrays
show_fig       = True   # applies to figures with statistical results
save_fig       = True   # applies to figures with statistical results
save_summed    = False  # applies to the contrasts images
save_averaged  = False  # averaged beta arrays
average_voxels = False  # CONTRASTS: If True, one value (array) across VOXELS.
                        # If both are False, the extracted roi array has shape (n_runs, n_voxels)
average_runs   = True   # If True, collapse runs and return a mean across runs.                  
remove_empty   = False  # Remove or not empty arrays (e.g.: If we do not average across voxels, 
					    # do we, when averaging across runs, include voxels that have zero 
					    # beta values or not?)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04. Specify experiment design and MRI info
subIDs = [5]
subID  = 5
anatID = 2
space  = "T1w" #TODO: What is the difference between T1w and T1wFOV?
task   = "timDev"
sesIDs = [2] # 2, 3, 4, 5, 6, 7
sessions = 2 # appears in the filenames 234567
acqIDs = ["BLOCK1"] # "FUNLOC", "BLOCK2", "BLOCK3", "BLOCK4"
blocks = "1" # appears in the filenames 1234

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 05. Specify project directories
# Core paths
# homePath  = Path("/home/mutrosa/mutrosa/Documents/projects/devLoc") # Citrix
homePath  = Path("/home/mutrosa/Documents/projects/devLoc")           # Local
baseDir, workDir, outDir = get_base_dirs(homePath, develop_mode, jobName, denoising)

# Define derived paths
dataPath   = outDir / "1stLevel"
out_2nd    = outDir / "2ndLevel"
out_1st    = dataPath / "visualization"
spmt_out   = baseDir / "visualization"
atlasPath  = homePath / "templates" / "resampled"
physioPath = homePath / "data_physio" / "raw"

# Preproc and filtered data paths
mriPath  = homePath / "data_MRI" / "derivatives" / f"NORDIC-{denoising}" / "derivatives" # path to preproc outputs
anatPath = mriPath / f"sub-{subID:02d}" / f"ses-{anatID:02d}" / "anat"
funcPath = mriPath / f"sub-{subID:02d}"
artPath  = homePath / "data_physio" / "artifacts" / f"NORDIC-{denoising}"
freesurfer_dir = mriPath / "sourcedata" / "freesurfer" / f"sub-{subID:02d}_ses-{anatID:02d}" / "mri"

# Create missing directories
for p in [outDir, out_2nd, out_1st, spmt_out, workDir]:
    p.mkdir(parents=True, exist_ok=True)

# Filenames and folders of the 1st level analysis output
# The 1st level results have to be resampled prior to visualization
con_name = "timDev-freqDev"            # contrast label
# TODO: boldref or T1w or T1wFOV (impacts resampling before visualization!) :below:
con_filename   = f"con_space-{space}_"   # image
beta_filename  = f"beta_space-{space}_"  # image
spmT_filename  = f"spmT_space-{space}_"  # image

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 06. Specify atlas and roi info
with open("rois.yaml", "r") as f:
    data = yaml.safe_load(f)
rois_cortical    = data["cortical"]
rois_subcortical = data["subcortical"]
rois = data["cortical"] | data["subcortical"]

# Get Sitek's subcortical atlas
atlas_subcor_name = f"sub-invivo_sub-{subID:02d}_ses-{anatID:02d}_space-{space}.nii.gz"
atlas_subcor_path = atlasPath / atlas_subcor_name

# Get FreeSurfer's parcellation: Destrieux Atlas
atlas_cor_name = f"aparc.a2009s+aseg_sub-{subID:02d}_ses-{anatID:02d}_NORDIC-{denoising}_space-{space}.nii.gz"
atlas_cor_path = atlasPath / atlas_cor_name

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