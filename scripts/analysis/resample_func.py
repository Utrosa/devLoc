#! /usr/bin/env python
# Time-stamp: <29-09-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI
# -----------------------------------------------------------------------------
'''
Resampling origin to boldref to T1 space with AFNI.
The anatomical reference - native T1 - is the same for all analysis outputs.
'''
from pathlib import Path
import grabber
import config as c
from objects import grab_objects
from utils import resample_img, compare_img

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BOLD functional images: from BOLDREF FOV to T1w FOV -------------------------
# Collect the files to resample and specify their resampled filename
for sesID in c.sesIDs:
    for acqID in c.acqIDs:

            # Grab files
            _, boldref_path, bold_path, _, _, _, _, _, T1w_path, \
            _, orig_to_boldref_path, boldref_to_T1w_path, _ = grab_objects(
                c.subID,  
                c.anatID, 
                c.homePath, 
                c.mriPath,
                c.artPath, 
                c.space,
                c.task, 
                ses=sesID,
                acq=acqID,
                run=None
            )

            # Select the bold and verify that it's only one
            if len(bold_path) == 1:
                bold = Path(bold_path[0])
            else:
                raise ValueError(f"\nFound too many/few bold files: {len(bold_path)}")

            # Select the bold reference and verify that it's only one
            if len(boldref_path) == 1:
                boldref = Path(boldref_path[0])
            else:
                raise ValueError(f"\nFound too many/few bold files: {len(boldref_path)}")

            # Select the transformation files and verify them
            if len(orig_to_boldref_path) == 1 and len(boldref_to_T1w_path) == 1:
                orig_to_boldref = Path(orig_to_boldref_path[0])
                boldref_to_T1w = Path(boldref_to_T1w_path[0])
            else:
                raise ValueError(f"\nFound too many/few transform files.\n"
                    f"\nOrig to boldref: {len(orig_to_boldref_path)}"
                    f"\nBoldref to T1: {len(boldref_to_T1w_path)}")

            # Change filenames to indicate change of space
            bold_boldref_name = bold.stem.replace(f"space-{c.space}", "space-orig") + bold.suffix
            bold_boldref_path = c.outDir / bold_boldref_name

            bold_T1_name = bold.stem.replace(f"space-{c.space}", c.resampled_stem) + bold.suffix
            bold_T1_path = c.outDir / bold_T1_name
            
            # Check that the resampled file does not already exist 
            if not bold_T1_path.exists():      
                if c.space == "T1w":
                    
                    # From origin to boldref
                    resample_img(
                        bold,
                        boldref, 
                        bold_boldref_path, 
                        "ants",
                        "NearestNeighbor",
                        orig_to_boldref)

                    # From boldref to T1w                    
                    # resample_img(
                    #     bold_boldref_path,
                    #     T1w_path,
                    #     bold_T1_path,
                    #     "ants",
                    #     "NearestNeighbor",
                    #     boldref_to_T1w)

                    # compare_img(bold, T1w_path, bold_new_path, verbose=c.verbose)

                else:
                    raise ValueError(f"What is the desired coordinate space? Not T1w nor MNI?")
            
            # Remove temporary files