#! /usr/bin/env python
# Time-stamp: <29-09-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI/nipypee
# Local:  conda activate localizer_fMRI/nipype
# -----------------------------------------------------------------------------
'''
Script written for looping over subjects. To be used after 1st level analysis.
Resample outputs from the 1st level analysis (beta images/t-values/contrasts)
from cropped restricted FoV space of functional scans to full space of the
anatomical image (T1w or MNI space).

Prerequisites: 
- install ANTs, nibabel, and nilearn
- run 1st level analysis
'''
import os
import tempfile
import subprocess
import config as c
import numpy as np
from nilearn import image
from nilearn.image import get_data
from nibabel import Nifti1Image, save
from utils import resample_img, compare_img

def clean_nans_infs(img, path):
    """
    Loads an image, replaces NaNs and Infs with 0 in memory, 
    and returns a new Nifti1Image object.
    """
    data = get_data(img)
    
    # Create a copy to avoid modifying the original data array if it's memory-mapped
    clean_data = data.copy() 
    clean_data[np.isnan(clean_data)] = 0
    clean_data[np.isinf(clean_data)] = 0

    # Create Nifti image
    clean_img = Nifti1Image(clean_data, img.affine, img.header)

    # Save a temporary file
    temp_fd, temp_path = tempfile.mkstemp(suffix = '.nii.gz', dir = str(path.parent))
    os.close(temp_fd)
    save(clean_img, temp_path)
    
    return temp_path

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 01. Choose the anatomical reference (it's the same for all analysis outputs).
if c.space == "T1w":
    template_img = c.anatPath / f"sub-{c.subID:02d}_ses-{c.anatID:02d}_desc-preproc_T1w.nii.gz"
else:
    template_img = c.homePath / "templates" / "tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz"
    print("\nAssumning analysis was done in MNI152NLin2009cAsym space.")

# 02. Collect the files to crop, specify their new filename and save to disk.
for sesID in c.sesIDs:
    for acqID in c.acqIDs:

            # All 1st-level results are in the same folder
            results_fold = c.dataPath / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}" / f"acq-{acqID}"
            
            # Check if the acquisition-specific folders exists
            if not results_fold.exists():
                results_fold = c.dataPath / f"sub-{c.subID:02d}" / f"ses-{sesID:02d}"
                print(f"\nFor ses-{sesID:02d} no subfolder 'acq-{acqID}' found.")

                if not results_fold.exists():
                    raise FileNotFoundError("The folder does not exist. Check the outputs of the analysis.")

            # Betas
            all_betas = list(results_fold.glob(f"beta_space-T1w*.nii*"))
            betas = [beta for beta in all_betas if c.resampled_stem not in beta.stem]
            for beta in betas:

                # Change filename to indicate resampling
                beta_new_name = beta.stem.replace(f"space-{c.space}", c.resampled_stem) + beta.suffix
                beta_new_path = beta.parent / beta_new_name
                
                # Replace NaNs and Infs with zeros
                beta_img = image.load_img(beta)                
                beta_temp_path = clean_nans_infs(beta_img, beta)
                del beta_img

                # Check that the resampled file does not already exist 
                if not beta_new_path.exists():      
                    if c.space == "T1w":
                        resample_img(beta_temp_path, template_img, beta_new_path, "nilearn", "nearest")
                        compare_img(beta_temp_path, template_img, beta_new_path, verbose=c.verbose)
                    else:
                        raise ValueError(f"What is the desired coordinate space? Not T1w nor MNI?")

                # Remove temporary file
                if os.path.exists(beta_temp_path):
                    os.remove(beta_temp_path)

            # Contrasts
            all_cons = list(results_fold.glob(f"con_space-T1w*.nii*"))
            cons = [con for con in all_cons if c.resampled_stem not in con.stem]
            for con in cons:

                # Change filename to indicate change to T1/MNI FOV space
                con_new_name = con.stem.replace(f"space-{c.space}", f"space-{c.space}FOV") + con.suffix
                con_new_path = con.parent / con_new_name
                
                # Replace NaNs and Infs with zeros
                con_img = image.load_img(con)                
                con_temp_path = clean_nans_infs(con_img, con)
                del con_img

                # Check that the resampled file does not already exist 
                if not con_new_path.exists():
                    if c.space == "T1w":     
                        resample_img(con_temp_path, template_img, con_new_path, "nilearn", "nearest")
                        compare_img(con_temp_path, template_img, con_new_path, verbose=c.verbose)
                    else:
                        raise ValueError(f"What is the desired coordinate space? Not T1w nor MNI?")
                
                # Remove temporary file
                if os.path.exists(con_temp_path):
                    os.remove(con_temp_path)
            
            # SPM t-images
            all_spmts = list(results_fold.glob(f"spmT_space-T1w*.nii*"))
            spmts = [spmt for spmt in all_spmts if c.resampled_stem not in spmt.stem]
            for spmt in spmts:

                # Change filename to indicate change to T1/MNI FOV space
                spmt_new_name = spmt.stem.replace(f"space-{c.space}", f"space-{c.space}FOV") + spmt.suffix
                spmt_new_path = spmt.parent / spmt_new_name
                
                # Replace NaNs and Infs with zeros
                spmt_img = image.load_img(spmt)                
                spmt_temp_path = clean_nans_infs(spmt_img, spmt)
                del spmt_img

                # Check that the resampled file does not already exist 
                if not spmt_new_path.exists():
                    if c.space == "T1w":
                        resample_img(spmt_temp_path, template_img, spmt_new_path, "nilearn", "nearest")
                        compare_img(spmt_temp_path, template_img, spmt_new_path, verbose=c.verbose)
                    else:
                        raise ValueError(f"What is the desired coordinate space? Not T1w nor MNI?")

                # Remove temporary file
                if os.path.exists(spmt_temp_path):
                    os.remove(spmt_temp_path)