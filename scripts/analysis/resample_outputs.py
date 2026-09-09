#! /usr/bin/env python
# Time-stamp: <07-09-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI
# -----------------------------------------------------------------------------
'''
Script written for looping over subjects. To be used after 1st level analysis.
A. Resample outputs from the 1st level analysis (beta images/t-values) from 
   restricted FoV space of functional scans to native space of the T1w image.

Prerequisites: 
- install ANTs, nibabel, and nilearn
- run 1st level analysis
'''
import subprocess
import config as c
from utils import resample_img, compare_img
# TODO: update so it does either conversion to MNI or T1, depending on space!
# TODO: the output from analysis should be which space (?!)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Beta, contrast, and spmT images: from BOLDREF FOV to T1w FOV ------------

# The anatomical reference is the same for all analysis outputs.
T1w = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_desc-preproc_T1w.nii.gz"

# Collect the files to resample and specify their resampled filename
for sesID in c.sesIDs:
    for acqID in c.acqIDs:

            # All 1st-level results are in the same folder
            results_fold = c.dataPath / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"

            # All results are subject to the same transform file
            from_boldref_to_T1w = c.funcPath / f"ses-{sesID:02d}" / "func" / f"sub-{c.subID:02d}_ses-{sesID:02d}_task-{c.task}_acq-{acqID}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
            if not from_boldref_to_T1w.exists():
                print(f"\nTransformation file is missing for sub-{c.subID}_ses-{sesID}_acq-{acqID}:\n {from_boldref_to_T1w}")
            
            # Betas
            betas = list(results_fold.glob(f"beta_space-{c.space}*.nii*"))
            for beta in betas:

                # Change filename to indicate change to T1/MNI FOV space
                beta_new_name = beta.stem.replace(f"space-{c.space}", f"space-{c.space}FOV") + beta.suffix
                beta_new_path = beta.parent / beta_new_name
                
                # Check that the resampled file does not already exist 
                if not beta_new_path.exists():      
                    resample_img(beta, T1w, beta_new_path, "ants", "NearestNeighbor", from_boldref_to_T1w)
                    compare_img(beta, T1w, beta_new_path)

            # Contrasts 
            cons = list(results_fold.glob(f"con_space-{c.space}*.nii*"))
            for con in cons:

                # Change filename to indicate change to T1/MNI FOV space
                con_new_name = con.stem.replace(f"space-{c.space}", f"space-{c.space}FOV") + con.suffix
                con_new_path = con.parent / con_new_name
                
                # Check that the resampled file does not already exist 
                if not con_new_path.exists():      
                    resample_img(con, T1w, con_new_path, "ants", "NearestNeighbor", from_boldref_to_T1w)
                    compare_img(con, T1w, con_new_path)
            
            # SPM t-images 
            spmts = list(results_fold.glob(f"spmt_space-{c.space}*.nii*"))
            for spmt in spmts:

                # Change filename to indicate change to T1/MNI FOV space
                spmt_new_name = spmt.stem.replace(f"space-{c.space}", f"space-{c.space}FOV") + spmt.suffix
                spmt_new_path = spmt.parent / spmt_new_name
                
                # Check that the resampled file does not already exist 
                if not spmt_new_path.exists():      
                    resample_img(spmt, T1w, spmt_new_path, "ants", "NearestNeighbor", from_boldref_to_T1w)
                    compare_img(spmt, T1w, spmt_new_path)