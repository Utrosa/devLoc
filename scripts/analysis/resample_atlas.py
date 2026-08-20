#! /usr/bin/env python
# Time-stamp: <21-06-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI
# -----------------------------------------------------------------------------
'''
A. Resample Sitek's in-vivo atlas to the resolution of the MNI template, used
   in preprocessing and data analyses, or to native space of the T1w image.

B. Resample images from 1st Level Analysis (results: beta images/t-values) from 
   restricted FoV space of functional scans to native space of the T1w image.

¡ DO NOT RUN THIS CODE AS A WHOLE but PER PART, DEPENDING ON YOUR GOAL !

Prerequisites: ANTs, nibabel, and nilearn.
'''
# -----------------------------------------------------------------------------
import subprocess
import nibabel as nib
from pathlib import Path
from nilearn.image import resample_to_img
from nibabel.orientations import aff2axcodes, io_orientation

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. PREPARATION -------------------------------------------------------------
# ------- Specify the subject and session identifiers
subID  = 5
anatID = 2 # Session in which the anatomical scan was collected
space   = "T1w"
denoising = True
sessions = [2, 3, 4, 5, 6, 7]
acqIDs = ["BLOCK1", "BLOCK2", "BLOCK3", "BLOCK4"] # "FUNCLOC"
task = "timDev" # "localizer"

# ------- Specify the project directories
homePath = Path("/home/mutrosa/mutrosa/Documents/devLoc")
mriPath  = homePath / "data_MRI" / "derivatives" / f"NORDIC-{denoising}" / "derivatives" / f"sub-{subID:02d}"
atlasPath = homePath / "templates"

# ------- Find and load the atlas
# Freesurfer (pre-covert to nifti with freesurfer) or Sitek (pre-resample to desired space)
# atlas_path =  atlasPath / f"aparc.a2009s+aseg_NORDIC-{denoising}_space-fsnative.nii.gz"
atlas_path = atlasPath / f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
atlas_img = nib.load(atlas_path)

# ------- Find and load the template images: MNI or T1w
# MNI_path   = homePath / "templates" / "tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz"
# MNI_img   = nib.load(MNI_path)   # MNI space in 1 mm
T1w_path   = mriPath / f"ses-{anatID:02d}" / "anat" / f"sub-{subID:02d}_ses-{anatID:02d}_desc-preproc_T1w.nii.gz"
T1w_img   = nib.load(T1w_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. RESAMPLING ATLAS to MNI 1x1x1 mmm ---------------------------------------
# Resample using nilearn (possible because the original atlas is in MNI space)
# resampled_atlas = resample_to_img(
#     atlas_img,
#     MNI_img,
#     interpolation='nearest'
# )
# Save
# out_path = out_dir / f"sub-invivo_resampled_to- MNI_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
# resampled_atlas.to_filename(out_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03. RESAMPLING ATLAS from MNI to T1 native space ----------------------------
# Apply a transform list by running an ANTs command ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -d 3: 3D images
# -i: input moving image (Atlas)
# -r: reference fixed image (T1w)
# -t: transform file (.h5 file for non-linear transformation and .)
# -n NearestNeighbor: critical for ROI labels to prevent interpolation artifacts
# -o: output path
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find the transform file
# transform_path = mriPath / f"ses-{anatID:02d}" / "anat" / f"sub-{subID:02d}_ses-{anatID:02d}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
#
# New name for the transformed file
# out_path = out_dir / f"sub-invivo_resampled_to-T1w_sub-{subID:02d}_ses-{anatID:02d}.nii.gz"
#
# Run the command
# cmd = [
#     "antsApplyTransforms",
#     "-d", "3",
#     "-i", str(atlas_path),
#     "-r", str(T1w_path),
#     "-t", str(transform_path),
#     "-n", "NearestNeighbor",
#     "-o", str(out_path)
# ]
# print(f"Running ANTs command: {' '.join(cmd)}")
# subprocess.run(cmd, check=True, capture_output=True, text=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04. RESAMPLING ATLAS from fsnative to T1 native space -----------------------
# Find the transform file
# transform_path = mriPath / f"ses-{anatID:02d}" / "anat" / f"sub-{subID:02d}_ses-{anatID:02d}_from-fsnative_to-T1w_mode-image_xfm.txt"
#
# New name for the transformed file
# out_path = atlasPath / f"aparc.a2009s+aseg_NORDIC-{denoising}_space-T1w.nii.gz"
#
# Run the command
# cmd = [
#     "antsApplyTransforms",
#     "-d", "3",
#     "-i", str(atlas_path),
#     "-r", str(T1w_path),
#     "-t", str(transform_path),
#     "-n", "NearestNeighbor",
#     "-o", str(out_path)
# ]
# print(f"Running ANTs command: {' '.join(cmd)}")
# subprocess.run(cmd, check=True, capture_output=True, text=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 05. RESAMPLING FUNC from BOLDREF FOV to T1w FOV -----------------------------
# for sesID in sessions:
#     for acqID in acqIDs:

#         # Functional image
#         func_path = mriPath / f"ses-{sesID:02d}" / "func" / f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_space-{space}_desc-preproc_bold.nii.gz"
#
#         # Find the transform file
#         transform_path = mriPath / f"ses-{sesID:02d}" / "func" / f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
#
#         # New name for the transformed file
#         out_path = mriPath / f"ses-{sesID:02d}" / "func" /f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_space-T1wFOV_desc-preproc_bold.nii.gz"

#         # Run the command
#         cmd = [
#             "antsApplyTransforms",
#             "-d", "3",
#             "-i", str(func_path),
#             "-r", str(T1w_path),
#             "-t", str(transform_path),
#             "-n", "NearestNeighbor",
#             "-o", str(out_path)
#         ]
#         print(f"\nRunning ANTs command: {' '.join(cmd)}")
#         subprocess.run(cmd, check=True, capture_output=True, text=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 06. RESAMPLING BETAS from BOLDREF FOV to T1w FOV ----------------------------
# conditions = list(range(1,12)); jobName = "when11where"
# conditions = [1, 2]; jobName = "whenwhat"
# beta_dir = homePath / "results" / jobName / f"NORDIC-{denoising}" / "1stLevel" / f"sub-{subID:02d}"
# for cond in conditions:
#     for sesID in sessions:
#         for acqID in acqIDs:

#             # Beta image from the analysis
#             beta_fold = beta_dir / f"ses-{sesID:02d}" / f"acq-{acqID}"
#             beta_name =  f"beta_space-boldref_{cond:04d}.nii"
#             beta_path = beta_fold / beta_name
#
#             # Find the transform file
#             transform_path = mriPath / f"ses-{sesID:02d}" / "func" / f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"            
#
#             # New name for the transformed file
#             out_path = beta_fold / f"beta_space-T1wFOV_{cond:04d}.nii"

#             # Run the command
#             cmd = [
#                 "antsApplyTransforms",
#                 "-d", "3",
#                 "-i", str(beta_path),
#                 "-r", str(T1w_path),
#                 "-t", str(transform_path),
#                 "-n", "NearestNeighbor",
#                 "-o", str(out_path)
#             ]
#             print(f"\nRunning ANTs command: {' '.join(cmd)}")
#             subprocess.run(cmd, check=True, capture_output=True, text=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 07. RESAMPLING con / spmT from BOLDREF FOV to T1w FOV -----------------------
res = "spmT" # con or spmT
jobName = "when11where" # whenwhat, when11where
res_dir = homePath / "results" / jobName / f"NORDIC-{denoising}" / "1stLevel" / f"sub-{subID:02d}"
for sesID in sessions:
    for acqID in acqIDs:

        # Contrast image from the analysis
        res_fold = res_dir / f"ses-{sesID:02d}" / f"acq-{acqID}"
        res_name =  f"{res}_space-boldref_0001.nii"
        res_path = res_fold / res_name

        # Find the transform file
        transform_path = mriPath / f"ses-{sesID:02d}" / "func" / f"sub-{subID:02d}_ses-{sesID:02d}_task-{task}_acq-{acqID}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"            
        
        # New name for the transformed file
        out_path = res_fold / f"{res}_space-T1wFOV_0001.nii"

        # Run the command
        cmd = [
            "antsApplyTransforms",
            "-d", "3",
            "-i", str(res_path),
            "-r", str(T1w_path),
            "-t", str(transform_path),
            "-n", "NearestNeighbor",
            "-o", str(out_path)
        ]
        print(f"\nRunning ANTs command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)

# 08. SANITY CHECKS --------------------------------------------------------------------------
# Check that the transformation has been correctly executed
resampled_img = nib.load(out_path)
resampled_shape  = resampled_img.shape
resampled_affine = resampled_img.affine

original_img = nib.load(atlas_path)
original_shape  = original_img.shape
original_affine = original_img.affine

template_img = T1w_img
template_shape  = template_img.shape
template_affine = template_img.affine

print(
    f"""Shape comparison:
- Original image shape   : {original_shape}
- Resampled image shape  : {resampled_shape}
- Template image shape   : {template_shape}
"""
)

print(
    f"""Affine comparison:
- Original image affine  : \n{original_affine}
- Resampled image affine : \n{resampled_affine}
- Template image affine  : \n{template_affine}
"""
)

print(
    f"""Axis direction codes comparison:
- Original image axcodes   : {aff2axcodes(original_affine)}
- Resampled image axcodes  : {aff2axcodes(resampled_affine)}
- Template image axcodes   : {aff2axcodes(template_affine)}
"""
)

print(
    f"""Input orientation comparison:
- Original image orientation  : \n{io_orientation(original_affine)}
- Resampled image orientation : \n{io_orientation(resampled_affine)}
- Template image orientation  : \n{io_orientation(template_affine)}
"""
)

print(
    f"""qform comparison:
- Original image qform        : {atlas_img.header.get_qform()[0]}
- Resampled image qform       : {resampled_img.header.get_qform()[0]}
- Template image qform        : {template_img.header.get_qform()[0]}
"""
)

print(
    f"""sform comparison:
- Original image sform        : {atlas_img.header.get_sform()[0]}
- Resampled image sform       : {resampled_img.header.get_sform()[0]}
- Template image sform        : {template_img.header.get_sform()[0]}
"""
)