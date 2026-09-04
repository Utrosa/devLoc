#! /usr/bin/env python
# Time-stamp: <31-08-2026 m.utrosa@bcbl.eu>
# Citrix: source activate localizer_fMRI
# Local:  conda activate localizer_fMRI

# Import python packages
import bids
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from nilearn.image import resample_to_img

# Import custom-made functions
import grabber

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

def extract_roi_array(subID, sesID, acqID, atlas, space, res_path, rois, out_dir, verbose, save, average_voxels):
    '''
    Extracts values from the specified regions of interest (ROIs).
    It's either one value per each voxel of the ROI or one value for the entire ROI.

    Parameters:
    - subID: integer number, identifying the participant. Only needded for the filename of files saved to disk.
    - sesID: integer number, identifying the session info. Only needded for the filename of files saved to disk.
    - acqID: string, identifying the functional MRI sequence. Only needded for the filename of files saved to disk.
    - atlas: string, path to an established atlas for auditory areas.
    - space: string, coordinate space of the input and output data (native T1w or MNI).
    - res_path: string, path to outputs of 1st Level Analysis.
    - rois: dictionary, specifying names, volume and atlas label of the ROIs.
    - out_dir: string, specifying the folder name for saving the results.
    - verbose: If True, prints affines and shape of atlas and result data in the terminal.
    - save: If True, saves the extracted roi arrays as a nifti files to disk.
    - average_voxels: If True, averages values across voxels of the ROI.

    Returns:
    - res_rois: dictionary, extracted value(s) per each ROI.
    - res_roi_paths: dictionary, paths to the extracted value(s) per each ROI.
    - res_affine: affine of the input nifti image (res_path).
    - Optionally: all extracted values are saved in out_dir.

    '''
    # Load atlas image
    atlas_img    = nib.load(atlas)
    atlas_data   = atlas_img.get_fdata()
    atlas_affine = atlas_img.affine

    # Load image from the analysis: beta image
    res_img    = nib.load(res_path)
    res_data   = res_img.get_fdata() 
    res_affine = res_img.affine
    
    # Compare shape, affine, qform, and sform
    if verbose:

        # Shape
        print("\natlas shape\n", atlas_img.shape)
        print("\ninput shape\n", res_img.shape)

        # Affines
        print("\natlas affine\n", atlas_affine)
        print("\ninput affine\n", res_affine)

        # Q form
        print("\n\nqform res\n", res_img.header.get_qform()[0])
        print("\nqform atlas\n",  atlas_img.header.get_qform()[0])

        # S form
        print("\n\nsform res\n", res_img.header.get_sform()[0])
        print("\nsform atlas\n",  atlas_img.header.get_sform()[0])

    # Extract values per ROIs
    res_rois = {}
    res_roi_paths = {}
    for name, roi in rois.items():

        # From the atlas, extract an array representing the roi (the mask)
        mask_data  = (atlas_data == roi['label']).astype(float) # 3D array with booleans -> floats

        # Count the voxels in the mask
        mask_size = np.sum(mask_data)

        # Apply the mask to find the ROI in result data
        res_array = res_data[mask_data > 0] # 1D array

        # Print the number of extracted voxels for the current ROI 
        if verbose:
            print(f"\nThere are {mask_size} voxels in the atlas for {name} region for result image:\n {res_path}")

        # Optionally save result as a zipped nifti file
        if save:
            res_masked = mask_data * res_data
            result_filename = f"sub-{subID:02d}_ses-{sesID:02d}_acq-{acqID}_roi-{name}_space-{space}.nii.gz"
            result_path = out_dir / result_filename
            res_roi_paths[name] = result_path
            nib.save(nib.Nifti1Image(res_masked, res_affine), result_path)

        # Append the extracted values for further analysis or visualization
        if average_voxels: # Collapse voxels: returns an array of shape (n_runs,)
            res_rois[name] = np.mean(res_array, axis=0)

        else: # Keep full data: (n_runs, n_voxels)
            res_rois[name] = res_array

    return res_rois, res_roi_paths, res_affine

def plot_violins(mask_paths, subID, sesID, acqIDs, out_dir, space, scale):

    rows = []
    for acq_name, roi_masks in mask_paths.items():
        
        for roi_name, roi_path in roi_masks.items():
            mask_img = nib.load(roi_path)
            vals = mask_img.get_fdata().flatten()
            vals = vals[vals != 0]
            if len(vals) < 5:
                print(f"Warning: very few voxels for {roi_name}, {acq_name}")
            rows.extend([{"ROI": roi_name, "acqID": acq_name, "values": v} for v in vals])
        
    df = pd.DataFrame(rows)
    print(df.head())

    for roi, group in df.groupby("ROI"):
        n_acq = len(acqIDs)
        fig, axes = plt.subplots(1, n_acq, figsize = (1.5 * n_acq, 8), sharey = True)
        
        if n_acq == 1:
            axes = [axes]
        
        for ax, acq in zip(axes, acqIDs):
            sub_df = group[group["acqID"] == acq]
            color_map = dict(zip(acqIDs, sns.color_palette("pastel", n_colors=len(acqIDs))))
            if not sub_df.empty:
                sns.violinplot(
                    y = "values",
                    data = sub_df,
                    ax = ax,
                    hue="ROI",
                    legend = False,
                    inner = "point",
                    cut = 0, 
                    palette = [color_map[acq]],
                    bw_adjust = 0.5
                )

                ax.set_title(f"{acq}", fontsize = 8)
                ax.set_xlabel("")
            ax.set_xticks([])
            if scale == True:
                ax.set_ylim(-5, 10)

        fig.suptitle(f"sub-{subID:02d}_ses-{sesID:02d}_roi-{roi}", fontsize = 12)
        fig.tight_layout()
        fig_name = f"sub-{subID:02d}_ses-{sesID:02d}_roi-{roi}_space-{space}_violins.png"
        fig_path = out_dir / fig_name
        plt.savefig(fig_path, dpi = 200, bbox_inches = "tight")
        plt.close(fig)

def plot_violins_average(betas, subID, sessions, blocks, plot_rois, n_cols, out_dir, space, scale, save, average_runs, average_voxels):
    """
    Plots the input data per ROI (subplots) and per condition (x axis categories).
    Depending on betas' structure, individual values in violin plots can be per run or per voxel.

    Parameters:
    - betas: nested dict with a list of beta arrays (single or multiple values) per ROI and condition.
    - subID: integer number, identifying the participant.
    - sessions: string of session numbers over which we are averaging.
    - blocks: string with blocks numbers over which we are averaging.
    - plot_rois: list of ROI labels. These will be the subplots of the figure.
    - n_cols: integer number subplots per column of the figure.
    - out_dir: string, specifying the folder name for saving the results as .nii.gz.
    - space: string, coordinate space of the input and output data (native T1w or MNI).
    - scale: If True, all subplots share the same y axis (scaled).
    - save: If True, saves the figure to disk.
    - average_runs: If True, "beta_array" is list with a single array of n_voxel values.
    - average_voxels: If True, "beta_array" is a list of n_runs integers.

    Returns:
    - Figure with violin subplots saved as .png
    """

    # Initialize a list to store values 
    rows = []

    # Iterate through each roi
    for roi in betas.keys():

        # Only plot data for the selected regions
        if roi in plot_rois:

            # Iterate through each condition
            for cond in betas[roi]:
                rows.extend([{"ROI": roi, "cond": cond, "values": v} for v in betas[roi][cond][0]])
        
        # Create a dataframe suitable for plotting      
        df = pd.DataFrame(rows)

    # Get unique values for coloring of violin plots (one per condition)
    conds = df["cond"].unique()
    violins = sns.color_palette("Set2", n_colors=len(conds))

    # Create a grid of subplots
    n_rows = int(np.ceil(len(plot_rois) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 10), sharey=True)
    axes = axes.flatten()

    # Plotting
    for i, roi in enumerate(plot_rois):
        ax = axes[i]
        
        # Filter data for the current ROI
        roi_data = df[df["ROI"] == roi]

        sns.violinplot(
            x = "cond",
            y = "values", 
            data = roi_data,
            ax = ax,
            hue = "cond",
            legend = False,
            inner = "point", # Show individual observations
            palette = violins,
            bw_adjust = 0.8,
            cut=0 # Limit the violin within the data range!
        )

        ax.set_title(f"{roi}", fontsize = 12)
        ax.set_xlabel("", fontsize = 12)
        ax.set_ylabel("", fontsize = 12)

        if scale == True:
            ax.set_ylim(-10, 10)

        fig.suptitle(f"sub-{subID:02d}", fontsize = 16, fontweight = "bold")
        fig.supxlabel("Timing Deviation [msec]", fontsize = 12)
        fig.supylabel("Beta Estimate [β]", fontsize = 12)
        fig.tight_layout()

    # Optionally save
    if save:
        fig_name = f"sub-{subID:02d}_ses-{sessions}_block-{blocks}_space-{space}_avgVox-{average_voxels}_avgRun-{average_runs}.png"
        fig_path = out_dir / fig_name
        plt.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Show the plot
    plt.show()

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

