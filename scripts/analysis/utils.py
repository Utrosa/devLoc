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

def check_for_missing_bunch_conditions(contrast_conditions, conditions_found, log_bunch, logfilepath):
    """
    The designs script includes functions that read events that occured during the experiment.
    It could happen that a designed event does not occur during a run.
    In that case, the designs script will return Bunch objects without that event (because the
    event is not in the log).

    The contrastor node in the GLM analysis with Nipype will fail if all runs do not have the
    same conditions. Therefore, we must add missing conditions to the bunch with this function.
    contrast_conditions: The levels of the manipulated variables in the experimental task that
                         are to be used in contrast estimation.
    conditions_found: The levels of the manipulated variables in the experimental task that
                      are present in logs of the run of the task.
    log_bunch: The bunch created by reading the logfiles.
    logfilepath: A Path to the log from which the log_bunch is made.
    """
    from nipype.interfaces.base import Bunch

    # Identify missing conditions and their positions
    missing = [con for ix, con in enumerate(contrast_conditions) if con not in conditions_found]

    if missing:
        logname = Path(logfilepath).name
        print(f"\nThe logfile with name {logname} is missing condition/s: {missing}.")

        # Build new lists that match the length of contrast_conditions
        new_onsets = [[] for _ in range(len(contrast_conditions))]
        new_durations = [[] for _ in range(len(contrast_conditions))]
        
        # Fill in existing data where conditions match
        # Where they don't match - leave empty lists
        for idx, cond in enumerate(conditions_found):
            if cond in contrast_conditions:
                target_idx = contrast_conditions.index(cond)
                
                # Find the corresponding onset and duration from the original bunch
                # Assuming original bunch has matching order for existing conditions
                try:
                    orig_idx = conditions_found.index(cond)
                    new_onsets[target_idx] = log_bunch.onsets[orig_idx]
                    new_durations[target_idx] = log_bunch.durations[orig_idx]
                except (IndexError, ValueError):
                    pass
        
        # Update the timfreq_bunch with the complete structure
        timfreq_bunch_corrected = Bunch(
            conditions=contrast_conditions,
            onsets=new_onsets,
            durations=new_durations)
    else:
        timfreq_bunch_corrected = log_bunch
        print("\nAll conditions are present, no changes needed.")

    return timfreq_bunch_corrected

def find_dev_group(delta_str, groups):
    """
    Determines the group name for a given delta string based on the groups dictionary.
    
    Args:
        delta_str (str): The delta value, potentially with a prefix (e.g., "p10", "n5", or "10").
        groups (dict): A dictionary where keys are group names and values are lists of integers.
        
    Returns:
        str|int: The name of the group if found, otherwise the parsed integer value.
    """
    for group_name, group_list in groups.items():
        if target_val in group_list:
            return group_name

    # If deviation not found in any group, return the raw value (or signed value)
    print(f"{target_val} not found in any group: {group_name}.")
    return target_val

def get_base_dirs(homePath, develop_mode, jobName, denoising):
    """Returns the core result directories based on mode."""
    base = homePath / ("results" if not develop_mode else "tests")
    work = base / f"work-{jobName}" / f"NORDIC-{denoising}"
    out  = base / jobName / f"NORDIC-{denoising}"
    return base, work, out

def compare_img(original_img_path, template_img_path, resampled_img_path, verbose=False):
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

    if verbose:
        print(
            f"""Shape comparison:
        - Original image shape   : {original_shape}
        - Resampled image shape  : {resampled_shape}
        - Template image shape   : {template_shape}
        """
        )

        print(
            f"""Affine comparison:
        - Original image affine  : \n{original_affine}\n
        - Resampled image affine : \n{resampled_affine}\n
        - Template image affine  : \n{template_affine}\n
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

def extract_roi_array(atlas, space, res_path, rois, out_dir, verbose, save, average_voxels):
    '''
    Extracts values from the specified regions of interest (ROIs).
    It's either one value per each voxel of the ROI or one value for the entire ROI.

    Parameters:
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
        res_masked = mask_data * res_data
        result_path = res_path.parent / f"roi-{name}_{res_path.stem}.nii.gz"
        res_roi_paths[name] = result_path
        if save:
            nib.save(nib.Nifti1Image(res_masked, res_affine), result_path)

        # Append the extracted values for further analysis or visualization
        if average_voxels: # Collapse voxels: returns an array of shape (n_runs,)
            res_rois[name] = np.mean(res_array, axis=0)

        else: # Keep full data: (n_runs, n_voxels)
            res_rois[name] = res_array

    return res_rois, res_roi_paths, res_affine

def plot_violins(mask_paths, subID, sesID, acqIDs, out_dir, space, scale):

    rows = []
    for acq_name, acq_masks in mask_paths.items():
        for ses_name, roi_masks in acq_masks.items():
            for roi_name, roi_path in roi_masks.items():
                mask_img = nib.load(roi_path)
                vals = mask_img.get_fdata().flatten()
                vals = vals[vals != 0]
                if len(vals) < 5:
                    print(f"Warning: very few voxels for {roi_name}, {acq_name}")
                rows.extend([{"ROI": roi_name, "acqID": acq_name, "sesID": ses_name, "values": v} for v in vals])
            
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

def plot_violins_average(betas, subID, plot_rois, plot_conf, out_dir, space, scale, save, show, average_runs, average_voxels):
    """
    Plots the input data per ROI (subplots) and per condition (x axis categories).
    Depending on betas' structure, individual values in violin plots can be per run or per voxel.

    Parameters:
    - betas: nested dict with a list of beta arrays (single or multiple values) per ROI and condition.
    - subID: integer number, identifying the participant.
    - plot_rois: list of ROI labels. These will be the subplots of the figure.
    - plot_conf: dictionary with plot configuration (the number of columns, ...).
    - out_dir: string, specifying the folder name for saving the results as .nii.gz.
    - space: string, coordinate space of the input and output data (native T1w or MNI).
    - scale: If True, all subplots share the same y axis (scaled).
    - save: If True, saves the figure to disk.
    - show: If True, shows the figure from terminal.
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
                data = betas[roi][cond][0]
                if np.ndim(data) == 0:
                    values_to_iterate = [data]
                else:
                    values_to_iterate = data
                rows.extend([{"ROI": roi, "cond": cond, "values": v} for v in values_to_iterate])
        
        # Create a dataframe suitable for plotting      
        df = pd.DataFrame(rows)

    # Get unique values for coloring of violin plots (one per condition)
    conds = df["cond"].unique()
    violins = sns.color_palette("Set2", n_colors=len(conds))

    # Create a grid of subplots
    n_rows = int(np.ceil(len(plot_rois) / plot_conf["cols"]))
    fig, axes = plt.subplots(n_rows, plot_conf["cols"], figsize=plot_conf["figsize"], sharey=True)
    axes = axes.flatten()

    # Plotting
    for i, roi in enumerate(plot_rois):
        ax = axes[i]
        
        # Filter data for the current ROI
        roi_data = df[df["ROI"] == roi]

        # Create violin shape
        sns.violinplot(
            data=roi_data,
            x="cond",
            y="values",
            hue="cond",
            palette=violins,
            inner="point", # show individual observations
            legend=False,
            cut=0, # limit the violin within the data range
            ax=ax
        )

        ax.set_title(f"{roi}", fontsize=plot_conf["subplot_fontsize"])
        ax.set_xlabel("", fontsize=plot_conf["subplot_fontsize"])
        ax.set_ylabel("", fontsize=plot_conf["subplot_fontsize"])

        # Add reference line at y=0
        ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8, zorder=1)

        # Optionally, limit the y axis
        if scale == True:
            ax.set_ylim(-10, 10)

        # Set figure title and shared axis labels
        fig.suptitle(
            f"sub-{subID:02d}, avgVox: {average_voxels}, avgRun: {average_runs}",
            fontsize=plot_conf["fig_fontsize"],
            fontweight="bold")
        fig.supxlabel(
            "Timing Deviation [msec]",
            fontsize=plot_conf["fig_fontsize"],
            fontweight="bold")
        fig.supylabel(
            "Beta Estimate [β]",
            fontsize=plot_conf["fig_fontsize"],
            fontweight="bold")
        fig.tight_layout()

    # Optionally save
    if save:
        fig_name = f"sub-{subID:02d}_space-{space}_avgVox-{average_voxels}_avgRun-{average_runs}.png"
        fig_path = out_dir / fig_name
        plt.savefig(fig_path, dpi=plot_conf["dpi"], bbox_inches="tight")

    # Show the plot
    if show:
        plt.show()
    else:
        plt.close(fig)

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
            interpolation,
            copy_header=True, # Copy the header of the input image to output image
            force_resample=True
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

def add_nuisance(bunch_dict, confounds_list, confounds_names):
    '''
    Parameters:
        bunch_dict: a list with Bunch objects, created by parsing logfiles of the experimental task.
        confounds_list: list of paths to filtered confounds (physiological regressors files - TAPAS!).
        confounds_names: column names of these confounds/regressors.

    Returns:
        List: a list of lists with bunch objects that include the regressors.
    '''
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from nipype.interfaces.base import Bunch

    design_bunch_list = []
    for bunch_log, conf_path in zip(bunch_dict, confounds_list):

        # Verify that the selected bunch and regressors refer to the same
        if not bunch_log.removesuffix("_events") == Path(conf_path).stem.removesuffix("_regressors"):
            raise ValueError(
                "The bunch object and regressors file do not correspond to the same sub/ses/acq:\n"
                f"Bunch: {bunch_log.removesuffix('_events')}\n"
                f"Regressors: {Path(conf_path).stem.removesuffix('_regressors')}"
            )

        # Read the confounds file
        all_confounds = pd.read_csv(
            conf_path,
            sep='\t',
            header=None,
            names=confounds_names,
            index_col=False
        )

        # Remove the last row of the dataframe which correspnds to the noise scan volume
        # This is not removed for physiological data!
        all_confounds = all_confounds.iloc[:-1]

        # Convert to the required format for SPM Bunch
        regressors = [all_confounds[col].tolist() for col in all_confounds.columns]
        
        # Regressor validation
        empty_or_nan_cols = []
        for i, col in enumerate(all_confounds.columns):
            series = all_confounds[col]
            
            # Check 1: Is the series empty?
            if series.empty:
                empty_or_nan_cols.append(col)
                continue
                
            # Check 2: Is the series full of NaNs?
            if series.isna().all():
                empty_or_nan_cols.append(col)

        if empty_or_nan_cols:
            print(f"Warning: The following regressors are empty or contain only NaNs and will be skipped: {empty_or_nan_cols}")
            
            # Filter them out
            valid_regressors = [r for i, r in enumerate(regressors) if all_confounds.columns[i] not in empty_or_nan_cols]
        else:
            print("All regressors contain valid data.")
            valid_regressors = regressors
        
        # Select the correct bunch
        design_bunch = bunch_dict[bunch_log]

        # Add regressors
        # https://nipype.readthedocs.io/en/1.11.0/api/generated/nipype.algorithms.modelgen.html
        # CHECK: They say they want "regressors" for the Bunch but with concatenate=False, modeler fails with that name
        # and it requires "regress"!
        design_bunch.regress = valid_regressors
        design_bunch.regressor_names = confounds_names

        # Append to list
        design_bunch_list.append(design_bunch)
    

    return design_bunch_list