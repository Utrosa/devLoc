#! /usr/bin/env python
# Time-stamp: <24-06-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI
import subprocess
import nibabel as nib
from pathlib import Path
from nilearn.image import resample_to_img

def resample_img(target, reference, output, method, interpolation, transform=""):
    """
    Resamples a target NifTi images to the space and resolution of the reference.
    The nilearn method can be used for resampling when:
        a.) both images are in the same coordinate space (MNI, T1w, or fsnative), and 
        b.) you only want to change the resolution (voxel size).

    The ANTs method has to be used when resampling to a different space and resolution.
    Here, you must select the correct transform file per subject/session/run:
        a.) from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
        b.) from-fsnative_to-T1w_mode-image_xfm.txt
        c.) from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt
    """
    input_path     = Path(target)
    reference_path = Path(reference)
    transform_path = Path(transform)
    output_path    = Path(output)

    if method == "nilearn":
    
        # Nilearn command
        input_img     = nib.load(input_path)
        reference_img = nib.load(reference_path)

        resampled_atlas = resample_to_img(
            input_img,
            reference_img,
            interpolation
        )
        resampled_atlas.to_filename(output_path)

    elif method == "ants":

        # ANTs command
        # -d 3: 3D images
        # -i: input moving image (Atlas)
        # -r: reference fixed image (T1w)
        # -t: transform file (.h5 file for non-linear transformation and .)
        # -n NearestNeighbor: critical for ROI labels to prevent interpolation artifacts
        # -o: output path
        cmd = [
            "antsApplyTransforms",
            "-d", "3",
            "-i", str(input_path),
            "-r", str(reference_path),
            "-t", str(transform_path),
            "-n", interpolation,
            "-o", str(output_path)
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"Erreur lors du resampling : {e.stderr}")
            return False
        except FileNotFoundError:
            print("Erreur: antsApplyTransforms n'est pas trouvé dans le PATH.")
            return False

def compare_img(original_img_path, template_img_path, resampled_img_path):
    """
    Compares three nifti images in shape, affines, and form.
    Useful for checking that resamping has been correctly executed
    """
    from nibabel.orientations import aff2axcodes, io_orientation

    original_img    = nib.load(original_img_path)
    original_shape  = original_img.shape
    original_affine = original_img.affine

    resampled_img    = nib.load(resampled_img_path)
    resampled_shape  = resampled_img.shape
    resampled_affine = resampled_img.affine

    template_img    = nib.load(template_img_path)
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
    - Original image qform        : {original_img.header.get_qform()[0]}
    - Resampled image qform       : {resampled_img.header.get_qform()[0]}
    - Template image qform        : {template_img.header.get_qform()[0]}
    """
    )

    print(
        f"""sform comparison:
    - Original image sform        : {original_img.header.get_sform()[0]}
    - Resampled image sform       : {resampled_img.header.get_sform()[0]}
    - Template image sform        : {template_img.header.get_sform()[0]}
    """
    )