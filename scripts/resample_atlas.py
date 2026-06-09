#! /usr/bin/env python
# Time-stamp: <12-05-2026 m.utrosa@bcbl.eu>
# conda activate localizer_fMRI
# -----------------------------------------------------------------------------
'''
Resample Sitek's in-vivo atlas to the resolution of the MNI template, used
in preprocessing and data analyses, or to native space of the T1w image, or 
restricted space of beta images.

DO NOT RUN THIS CODE AS A WHOLE but PER PART, DEPENDING ON YOUR GOAL.

Prerequisites: ANTs, nibabel, and nilearn
DOI: 10.7554/eLife.48932
'''
# -----------------------------------------------------------------------------
import subprocess
import nibabel as nib
from pathlib import Path
from nilearn.image import resample_to_img
from nibabel.orientations import aff2axcodes, io_orientation

# 01. PREPARATION -------------------------------------------------------------
# Specify the subject and session identifiers.
subID = 5
sesID = 2
space = "T1w"
denoising = True

# Specify the project directory.
homePath   = Path("/home/mutrosa/mutrosa/Documents/devLoc")

# Find the atlas, MNI, T1 and transform files.
atlas_path = homePath / "templates" / f"sub-invivo_resampled_to-{space}_sub-{subID:02d}_ses-{sesID:02d}.nii.gz"

MNI_path   = homePath / "templates" / "tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz"
T1w_path   = homePath / "data_MRI" / "derivatives" / f"NORDIC-{denoising}" / "derivatives" / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / "anat" / f"sub-{subID:02d}_ses-{sesID:02d}_desc-preproc_T1w.nii.gz"
transform_path = homePath / "data_MRI" / "derivatives" / f"NORDIC-{denoising}" / "derivatives" / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / "anat" / f"sub-{subID:02d}_ses-{sesID:02d}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"

beta_path  = homePath / "results" / "timDev_abs" / f"NORDIC-{denoising}" / "1stLevel" / f"sub-{subID:02d}" / f"ses-{sesID:02d}" / "acq-BLOCK1" / f"beta_space-{space}_0011.nii"# can be any beta in the same space as atlas

# Specify where to save the resampled atlases
out_dir = homePath / "templates"
out_dir.mkdir(exist_ok=True)

# Load images
atlas_img = nib.load(atlas_path) 
MNI_img   = nib.load(MNI_path)   # MNI space in 1 mm
T1w_img   = nib.load(T1w_path)   # Native space
beta_img  = nib.load(beta_path)  # Native restricted space

# 02. RESAMPLING to MNI 1x1x1 mmm -------------------------------------------------------------
# Resample using nilearn (possible because atlas is in MNI space)
resampled_atlas = resample_to_img(
    atlas_img,
    MNI_img,
    interpolation='nearest'
)

# Save
out_path = out_dir / f"sub-invivo_resampled_to- MNI_sub-{subID:02d}_ses-{sesID:02d}.nii.gz"
resampled_atlas.to_filename(out_path)

# 03. RESAMPLING to T1 native space ------------------------------------------------------------
# Apply a transform list by running an ANTs command
# -d 3: 3D images
# -i: input moving image (Atlas)
# -r: reference fixed image (T1w)
# -t: transform file (.h5 file for non-linear transformation and .)
# -n NearestNeighbor: critical for ROI labels to prevent interpolation artifacts
# -o: output path
out_path = out_dir / f"sub-invivo_resampled_to-T1w_sub-{subID:02d}_ses-{sesID:02d}.nii.gz"
cmd = [
    "antsApplyTransforms",
    "-d", "3",
    "-i", str(atlas_path),
    "-r", str(T1w_path),
    "-t", str(transform_path),
    "-n", "NearestNeighbor",
    "-o", str(out_path)
]
print(f"Running ANTs command: {' '.join(cmd)}")

try:
    # Run the command
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    print(f"Success! Resampled atlas saved to: {out_path}")
    
    # Verify the output exists and has the correct shape
    if out_path.exists():
        img = nib.load(out_path)
        t1w = nib.load(T1w_path)
        print(f"Output shape: {img.shape} (Matches T1w: {t1w.shape})")
    else:
        print("Error: Output file was not created.")

except subprocess.CalledProcessError as e:
    print(f"Error running antsApplyTransforms: {e}")
    print(f"stderr: {e.stderr}")
except FileNotFoundError:
    print("Error: 'antsApplyTransforms' command not found. Ensure your conda environment is activated.")

# 04. RESAMPLING to BETA ----------------------------------------------------------------------
# Resample using nilearn (possible because atlas is in MNI space)
resampled_atlas = resample_to_img(
    atlas_img, # Atlas resampled to native space
    beta_img,
    interpolation='nearest'
)

# Save
out_path = out_dir / f"sub-invivo_resampled_to-T1w_desc-betaNORDIC{denoising}_sub-{subID:02d}_ses-{sesID:02d}.nii.gz"
resampled_atlas.to_filename(out_path)

# 05. SANITY CHECKS --------------------------------------------------------------------------
# Check that the transformation has been correctly executed
resampled_atlas_img = nib.load(out_path)

original_shape  = atlas_img.shape
original_affine = atlas_img.affine

resampled_shape  = resampled_atlas_img.shape
resampled_affine = resampled_atlas_img.affine

MNI_shape  = MNI_img.shape
MNI_affine = MNI_img.affine

T1w_shape  = T1w_img.shape
T1w_affine = T1w_img.affine

print(
    f"""Shape comparison:
- Original atlas image shape  : {original_shape}
- Resampled atlas image shape : {resampled_shape}
- MNI template image shape    : {MNI_shape}
- T1w image shape             : {T1w_shape}
"""
)

print(
    f"""Affine comparison:
- Original atlas image affine  : \n{original_affine}
- Resampled atlas image affine : \n{resampled_affine}
- MNI template image affine    : \n{MNI_affine}
- T1w image affine             : \n{T1w_affine}
"""
)

print(
    f"""Axis direction codes comparison:
- Original atlas image axcodes  : {aff2axcodes(original_affine)}
- Resampled atlas image axcodes : {aff2axcodes(resampled_affine)}
- MNI template image axcodes    : {aff2axcodes(MNI_affine)}
- T1w image axcodes             : {aff2axcodes(T1w_affine)}
"""
)

print(
    f"""Input orientation comparison:
- Original atlas image orientation  : \n{io_orientation(original_affine)}
- Resampled atlas image orientation : \n{io_orientation(resampled_affine)}
- MNI template image orientation    : \n{io_orientation(MNI_affine)}
- T1w image orientation             : \n{io_orientation(T1w_affine)}
"""
)

print(
    f"""qform comparison:
- Original atlas image qform  : {atlas_img.header.get_qform()[0]}
- Resampled atlas image qform : {resampled_atlas_img.header.get_qform()[0]}
- MNI template image qform    : {MNI_img.header.get_qform()[0]}
- T1w image qform             : {T1w_img.header.get_qform()[0]}
"""
)

print(
    f"""sform comparison:
- Original atlas image sform  : {atlas_img.header.get_sform()[0]}
- Resampled atlas image sform : {resampled_atlas_img.header.get_sform()[0]}
- MNI template image sform    : {MNI_img.header.get_sform()[0]}
- T1w image sform             : {T1w_img.header.get_sform()[0]}
"""
)