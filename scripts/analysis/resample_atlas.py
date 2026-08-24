#! /usr/bin/env python
# Time-stamp: <24-06-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI
# -----------------------------------------------------------------------------
'''
Script written for looping over subjects.
A. Resample Sitek's in-vivo atlas to the resolution of the MNI template used
   in preprocessing and data analyses, and then to T1w native space.

B. Resample Freesurfer's reconall atlas to T1w native space.

C. Resample outputs from the 1st level analysis (beta images/t-values) from 
   restricted FoV space of functional scans to native space of the T1w image.

Prerequisites: 
- install ANTs, nibabel, and nilearn
- download Sitek's atlas and MNI template
'''
import config as c
import subprocess
from resample_utils import resample_img, compare_img

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 00. PREPARATION -------------------------------------------------------------
# ------- Specify the atlas-specific project directories
tempPath = c.homePath / "templates"
outAtlas = tempPath / "resampled"
outAtlas.mkdir(parents=True, exist_ok=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Sitek's in-vivo atlas: from MNI to T1 native space ----------------------
# Sitek's original atlas has resolution 0.5 mm isotropic (2009b!)
# MNI152NLAsym 0.5mm (b) -> MNI152NLAsym 1mm (c)
sitek_05mm = tempPath / f"sub-invivo_MNI_rois.nii.gz"
MNI_1mm    = tempPath / "tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz"
sitek_1mm  = outAtlas / f"sub-invivo_MNI152NLin2009cAsym_res-01.nii.gz"
if not sitek_1mm.exists():
    resample_img(sitek_05mm, MNI_1mm, sitek_1mm, "nilearn", "nearest")
    compare_img(sitek_05mm, MNI_1mm, sitek_1mm)

# MNI152NLAsym 1mm (c) -> T1
T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_desc-preproc_T1w.nii.gz"
from_MNI152NLin2009cAsym_to_T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
sitek_T1w = outAtlas / f"sub-invivo_T1w_sub-{c.subID:02d}_ses-{c.anatID:02d}.nii.gz"
if not sitek_T1w.exists():
    resample_img(sitek_1mm, T1w, sitek_T1w, "ants", "NearestNeighbor", from_MNI152NLin2009cAsym_to_T1w)
    compare_img(sitek_1mm, T1w, sitek_T1w)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Freesurfer's atlas: from fsnative to T1 native space --------------------
freesurfer_orig = tempPath / f"aparc.a2009s+aseg_sub-{c.subID:02d}_ses-{c.anatID:02d}_NORDIC-{c.denoising}_space-fsnative.nii.gz"
if not freesurfer_orig.exists():
    # Transform to nifti with mri_convert command  
    freesurfer_mgz = c.freesurfer_dir / "aparc.a2009s+aseg.mgz"
    if not freesurfer_mgz.exists():
        raise FileNotFoundError(f"Freesurfer source file not found: {freesurfer_mgz}")
    cmd = ["mri_convert", str(freesurfer_mgz), str(freesurfer_orig)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    
from_fsnative_to_T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_from-fsnative_to-T1w_mode-image_xfm.txt"
freesurfer_T1w = outAtlas / f"aparc.a2009s+aseg_sub-{c.subID:02d}_ses-{c.anatID:02d}_NORDIC-{c.denoising}_space-T1w.nii.gz"
if not freesurfer_T1w.exists():
    resample_img(freesurfer_orig, T1w, freesurfer_T1w, "ants", "NearestNeighbor", from_fsnative_to_T1w)
    compare_img(freesurfer_orig, T1w, freesurfer_T1w)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 03. Functional bold scans: from BOLDREF FOV to T1w FOV ----------------------
# TODO: Is this necessary? The transformed file show a white square when opening with fsleyes
for sesID in c.sesIDs:
    for acqID in c.acqIDs:
        func_orig = c.funcPath / f"ses-{sesID:02d}" / "func" / f"sub-{c.subID:02d}_ses-{sesID:02d}_task-{c.task}_acq-{acqID}_space-T1w_desc-preproc_bold.nii.gz"
        from_boldref_to_T1w = c.funcPath / f"ses-{sesID:02d}" / "func" / f"sub-{c.subID:02d}_ses-{sesID:02d}_task-{c.task}_acq-{acqID}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
        func_T1FOV = c.funcPath / f"ses-{sesID:02d}" / "func" / f"sub-{c.subID:02d}_ses-{sesID:02d}_task-{c.task}_acq-{acqID}_space-T1wFOV_desc-preproc_bold.nii.gz"

#         if not func_T1FOV.exists():
#             resample_img(func_orig, T1w, func_T1FOV, "ants", "NearestNeighbor", from_boldref_to_T1w)
#             compare_img(func_orig, T1w, func_T1FOV)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 04. Beta, contrast, and spmT images: from BOLDREF FOV to T1w FOV ------------
# TODO: Is this necessary? What does the shift achieve?
for cond in c.conditions_int:
    for sesID in c.sesIDs:
        for acqID in c.acqIDs:

            # All 1st-Level images are in the same folder and subject to the same transform
            # TODO: What is the correct resampling reference image here? 
            results_fold = c.data_1stLevel / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
            from_boldref_to_T1w = c.funcPath / f"ses-{sesID:02d}" / "func" / f"sub-{c.subID:02d}_ses-{sesID:02d}_task-{c.task}_acq-{acqID}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"

            # Betas
            beta_name =  f"beta_space-T1w_{cond:04d}.nii"
            beta_orig = results_fold / beta_name
            beta_T1FOV = results_fold / f"beta_space-T1wFOV_{cond:04d}.nii"

            # if not beta_T1FOV.exists():
                # resample_img(beta_orig, T1w, beta_T1FOV, "ants", "NearestNeighbor", from_boldref_to_T1w)
                # compare_img(beta_orig, T1w, beta_T1FOV)
                # subprocess.run(["freeview", str(beta_orig), str(T1w), str(beta_T1FOV)])
           
            # Contrasts 
            con_name =  "con_space-boldref_0001.nii"
            con_orig = results_fold / con_name
            con_T1FOV = results_fold / "con_space-T1wFOV_0001.nii"
            # TODO: add transform code 
            
            # SPM t-images 
            # TODO: check if the name is correct 
            spmt_name =  "spmt_space-boldref_0001.nii"
            spmt_orig = results_fold / spmt_name
            spmt_T1FOV = results_fold / "spmt_space-T1wFOV_0001.nii"
            # TODO: add transform code 