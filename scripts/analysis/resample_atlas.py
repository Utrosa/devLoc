#! /usr/bin/env python
# Time-stamp: <07-09-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI
# -----------------------------------------------------------------------------
'''
Script written for looping over subjects.
A. Resample Sitek's in-vivo atlas to the resolution of the MNI template used
   in preprocessing and data analyses, then to T1w native space (non-linear).

B. Resample Freesurfer's reconall atlas to T1w native space (linear), then to
   and MNI space (non-linear).

Prerequisites: 
- install ANTs, nibabel, and nilearn
- download Sitek's atlas and MNI template
'''
import shutil
import config as c
import subprocess
from utils import resample_img, compare_img

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
sitek_1mm  = outAtlas / f"sub-invivo_sub-{c.subID:02d}_ses-{c.anatID:02d}_space-MNI152NLin2009cAsym.nii.gz"
if not sitek_1mm.exists():
    resample_img(sitek_05mm, MNI_1mm, sitek_1mm, "nilearn", "nearest")
    if c.verbose:
        compare_img(sitek_05mm, MNI_1mm, sitek_1mm)

# MNI152NLAsym 1mm (c) -> T1
T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_desc-preproc_T1w.nii.gz"
from_MNI152NLin2009cAsym_to_T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
sitek_T1w = outAtlas / f"sub-invivo_sub-{c.subID:02d}_ses-{c.anatID:02d}_space-T1w.nii.gz"
if not sitek_T1w.exists():
    resample_img(sitek_1mm, T1w, sitek_T1w, "ants", "NearestNeighbor", from_MNI152NLin2009cAsym_to_T1w)
    if c.verbose:
        compare_img(sitek_1mm, T1w, sitek_T1w)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 02. Freesurfer's atlas: from fsnative to T1 native space --------------------
# Find the freesurfer atlas and move it to templates
# Note, recon-all will parcellate the individual subject’s brain according to 
#   the Desikan-Killiany atlas (aparc+aseg.mgz), and
#   the Destrieux atlas (aparc.a2009s+aseg.mgz; more parcellated).
freesurfer_nii = tempPath / f"aparc.a2009s+aseg_sub-{c.subID:02d}_ses-{c.anatID:02d}_NORDIC-{c.denoising}_space-fsnative.nii.gz"
freesurfer_mgz = c.freesurfer_dir / "aparc.a2009s+aseg.mgz"
freesurfer_tmp = tempPath / freesurfer_mgz.name
if not freesurfer_tmp.exists():
    if not freesurfer_mgz.exists():
        raise FileNotFoundError(f"Freesurfer source file not found: {freesurfer_mgz}")    
    shutil.copy2(freesurfer_mgz, freesurfer_tmp) # preserve file metadata

# Use freesurfer's utility function https://surfer.nmr.mgh.harvard.edu/fswiki/mri_convert)
# to covert the recon-all parcellation atlases from MGH to NifTi format.
if not freesurfer_nii.exists():
    cmd = ["mri_convert", str(freesurfer_tmp), str(freesurfer_nii)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)

# From fsnative to T1w
from_fsnative_to_T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_from-fsnative_to-T1w_mode-image_xfm.txt"
freesurfer_T1w = outAtlas / f"aparc.a2009s+aseg_sub-{c.subID:02d}_ses-{c.anatID:02d}_NORDIC-{c.denoising}_space-T1w.nii.gz"
if not freesurfer_T1w.exists():
    resample_img(freesurfer_nii, T1w, freesurfer_T1w, "ants", "NearestNeighbor", from_fsnative_to_T1w)
    if c.verbose:
        compare_img(freesurfer_nii, T1w, freesurfer_T1w)

# From T1w to MNI
from_T1w_to_MNI = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"
freesurfer_MNI  = outAtlas / f"aparc.a2009s+aseg_sub-{c.subID:02d}_ses-{c.anatID:02d}_NORDIC-{c.denoising}_space-MNI152NLin2009cAsym.nii.gz"
if not freesurfer_MNI.exists():
    resample_img(freesurfer_T1w, MNI_1mm, freesurfer_MNI, "ants", "NearestNeighbor", from_T1w_to_MNI)
    if c.verbose:
        compare_img(freesurfer_T1w, MNI_1mm, freesurfer_MNI)