#! /usr/bin/env python
# Time-stamp: <2026-15-06 m.utrosa@bcbl.eu>
'''
Grabs objects needed for other scripts.
Grabs functional and anatomical files in the specified space.
'''

def grab_objects(subID, sesID, anatID, homePath, mriPath, artPath, space, acqID, task, run):
	"""
	Locate functional and anatomical objects and returns a tuple of filepaths and
	TR based on the specified subject, session, acquisition, and run.

	Args:
		subID (int): The subject identifier.
		sesID (str): Session identifier.
		anatID (int): The identifier of the session in which the anatomical image was obtained.
		homePath (str): The base directory path.
		mriPath (str): The path to the MRI data.
		artPath (str): ThE path to the preprocessed physiological data.
		space (str): The target space for the files (e.g., MNI152NLin2009cAsym or T1w).
		acqID (str): Acquisition name.
		task (str): Name of the experimental task that the subject was doing.
		run (bool/str): If False, no run identifier is used to locate files. If str (True),
						that run identifier is used to find files.

	Returns:
		tuple: strings (filepaths) and TR (float).

	Note:
		Ensure that session, acquisition, and run match an existing scenario.
		Transformation files (xfm) are assumed to be returned in alphabetical order by the grabber.
	"""
	
	import bids
	import grabber
	import warnings
	from pathlib import Path

	if run:

		# If run is not False, then it is run identifier (a string).
		runID = run		

		# Main paths
		homePath = Path(homePath)
		mriPath  = Path(mriPath)
		artPath  = Path(artPath)

		# -------------- 01 Set up layouts -------------- 
		logpath   = homePath / "data_logs" / "bids"
		logLayout = bids.layout.BIDSLayout(logpath, validate=False)
		mriLayout = bids.layout.BIDSLayout(mriPath, validate=False, derivatives=True)
		artLayout = bids.layout.BIDSLayout(artPath, validate=False)
		
		# -------------- 02 Configuration -------------- 
		log_conf  = grabber.define_grabconf(subID, sesID, "events", "tsv", task = task, acquisition = acqID, run = runID)
		bold_conf = grabber.define_grabconf(subID, sesID, "bold", "nii.gz", task = task, acquisition = acqID, run = runID, space = space)
		mask_conf = grabber.define_grabconf(subID, sesID, "mask", "nii.gz", task = task, acquisition = acqID, run = runID, space = space)
		conf_conf  = grabber.define_grabconf(subID, sesID, "confounds", "tsv", acquisition = acqID, run = runID)	
		reg_conf = grabber.define_grabconf(subID, sesID, "regressors", "tsv", acquisition = acqID, run = runID)
		movpar_conf = grabber.define_grabconf(subID, sesID, "movpar", "txt", acquisition = acqID, run = runID)
		out_conf  = grabber.define_grabconf(subID, sesID, "outliers",  "txt", acquisition = acqID, run = runID)
		T1w_conf  = grabber.define_grabconf(subID, anatID, "T1w",  "nii.gz")
		T1w_to_MNI_conf = grabber.define_grabconf(subID, anatID, "xfm",  "h5")
		boldref_to_T1w_conf  = grabber.define_grabconf(subID, anatID, "xfm",  "txt")

		# -------------- 03 Grabbing files --------------
		log_object  = grabber.grab_BIDS_object(logpath, logLayout, log_conf)
		bold_object = grabber.grab_BIDS_object(mriPath, mriLayout, bold_conf)
		mask_object = grabber.grab_BIDS_object(mriPath, mriLayout, mask_conf)
		conf_object = grabber.grab_BIDS_object(artPath, artLayout, conf_conf) # selected confounds
		reg_object  = grabber.grab_BIDS_object(artPath, artLayout, reg_conf)  # only BIOPAC
		movpar_path = movpar_object[0].path # only the trans & rot parameters
		out_object    = grabber.grab_BIDS_object(artPath, artLayout, out_conf)  # motion outliers as detected by fMRIPrep
		T1w_object    = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_conf)
		T1w_to_MNI_object      = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_to_MNI_conf)
		orig_to_boldref_object = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)
		boldref_to_T1w_object  = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)

		# -------------- 04 Verification & Warnings --------------
		        
		# Check for missing files
		if len(log_object) == 0:
			raise ValueError(f"No log file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}, run-{runID:02d}.")
		if len(bold_object) == 0:
			raise ValueError(f"No bold file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}, run-{runID:02d}.")
		if len(mask_object) == 0:
			raise ValueError(f"No mask file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}, run-{runID:02d}.")
		if len(T1w_object) == 0:
			raise ValueError(f"No T1w file found for sub-{subID:02d}, ses-{anatID:02d}.")
		if len(conf_object) == 0:
			raise ValueError(f"No confounds file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}, run-{runID:02d}.")
		if len(reg_object) == 0:
			raise ValueError(f"No TAPAS regressors file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}, run-{runID:02d}.")
		if len(movpar_object) == 0:
			raise ValueError(f"No movement parameters file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")
		if len(out_object) == 0:
			raise ValueError(f"No outliers file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}, run-{runID:02d}.")
		if len(T1w_to_MNI_object) == 0:
			raise ValueError(f"No T1w to MNI transform file found sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}, run-{runID:02d}")
		if len(orig_to_boldref_object) == 0 or len(boldref_to_T1w_object) == 0:
			raise ValueError(f"No boldref to T1w transform files found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}, run-{runID:02d}")

		# Check for multiple files (ambiguity)
		if len(log_object) > 1:
			raise ValueError(
				f"Found more than one file:\n{log_object}."
				"\nCheck your LOG folder and logfile import steps."
				)
		
		if len(bold_object) > 1 or len(mask_object) > 1:
			raise ValueError(
				"Found more than one file in one of the following\n:"
				f"{bold_object}, or\n {mask_object}."
				"\nPlease check your MRI folder and MRI data import steps."
				)

		if len(conf_object) > 1 or len(movpar_object) > 1 or len(out_object) > 1:
			raise ValueError(
				"Found more than one file in one of the following\n:"
				f"{conf_object},\n {movpar_object}, or\n{out_object}."
				"Check your PHYSIO folder and physiological data import steps."
				)
	
		if len(reg_object) > 1:
			raise ValueError(
				"Found more than one file in one of the following\n:"
				f"{reg_object}."
				"Check your TAPAS preprocessing and regressors estimation steps."
				)

		# Warnings for T1 and transform files
		if len(T1w_object) > 1:
			warnings.warn(
				"Multiple anatomical files found. Taking the first one."
			)
		if len(orig_to_boldref_object) > 1 or len(boldref_to_T1w_object) > 1:
			warnings.warn(
                "Multiple transformation files found. Assuming alphabetical order: "
                "orig_to_boldref is the second file, boldref_to_T1w is the first. "
                f"Found: orig={orig_to_boldref_object}, boldref={boldref_to_T1w_object}"
            )
		
		# -------------- 05 Grabing filepaths and Updating --------------
		log_path    = log_object[0].path
		bold_path   = bold_object[0].path
		mask_path   = mask_object[0].path
		conf_path   = conf_object[0].path # selected confounds
		reg_path    = reg_object[0].path  # only BIOPAC
		movpar_path = movpar_object[0].path
		out_path    = out_object[0].path  # motion outliers as detected by fMRIPrep
		T1w_path    = T1w_object[0].path
		T1w_to_MNI_path      = T1w_to_MNI_object[1].path
		orig_to_boldref_path = orig_to_boldref_object[1].path
		boldref_to_T1w_path  = boldref_to_T1w_object[0].path

		# Extract repetition time with PyBIDS methods [sec]
		TR = bold_object[0].get_metadata()['RepetitionTime']
		
	else:
		# Main paths
		homePath = Path(homePath)
		mriPath  = Path(mriPath)
		artPath  = Path(artPath)

		# -------------- 01 Set up layouts -------------- 
		logpath   = homePath / "data_logs" / "bids"
		logLayout = bids.layout.BIDSLayout(logpath, validate=False)
		mriLayout = bids.layout.BIDSLayout(mriPath, validate=False, derivatives=True)
		artLayout = bids.layout.BIDSLayout(artPath, validate=False)
		
		# -------------- 02 Configuration -------------- 
		log_conf  = grabber.define_grabconf(subID, sesID, "events", "tsv", task = task, acquisition = acqID)
		bold_conf = grabber.define_grabconf(subID, sesID, "bold", "nii.gz", task = task, acquisition = acqID, space = space)
		mask_conf = grabber.define_grabconf(subID, sesID, "mask", "nii.gz", task = task, acquisition = acqID, space = space)
		conf_conf = grabber.define_grabconf(subID, sesID, "confounds", "tsv", acquisition = acqID)	
		reg_conf  = grabber.define_grabconf(subID, sesID, "regressors", "tsv", acquisition = acqID)
		movpar_conf = grabber.define_grabconf(subID, sesID, "movpar", "txt", acquisition = acqID)
		out_conf  = grabber.define_grabconf(subID, sesID, "outliers",  "txt", acquisition = acqID)
		T1w_conf  = grabber.define_grabconf(subID, anatID, "T1w",  "nii.gz")
		T1w_to_MNI_conf = grabber.define_grabconf(subID, anatID, "xfm",  "h5")
		boldref_to_T1w_conf  = grabber.define_grabconf(subID, anatID, "xfm",  "txt")

		# -------------- 03 Grabbing files --------------
		log_object   = grabber.grab_BIDS_object(logpath, logLayout, log_conf)
		bold_object  = grabber.grab_BIDS_object(mriPath, mriLayout, bold_conf)
		mask_object  = grabber.grab_BIDS_object(mriPath, mriLayout, mask_conf)
		conf_object  = grabber.grab_BIDS_object(artPath, artLayout, conf_conf)  # selected confounds
		reg_object   = grabber.grab_BIDS_object(artPath, artLayout, reg_conf)   # only BIOPAC
		movpar_object = grabber.grab_BIDS_object(artPath, artLayout, movpar_conf) # only the trans & rot parameters
		out_object    = grabber.grab_BIDS_object(artPath, artLayout, out_conf)    # motion outliers as detected by fMRIPrep
		T1w_object    = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_conf)
		T1w_to_MNI_object      = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_to_MNI_conf)
		orig_to_boldref_object = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)
		boldref_to_T1w_object  = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)

		# -------------- 04 Verification & Warnings --------------
		# Check for missing files
		if len(log_object) == 0:
			raise ValueError(f"No log file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")
		if len(bold_object) == 0:
			raise ValueError(f"No bold file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}.")
		if len(mask_object) == 0:
			raise ValueError(f"No mask file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}.")
		if len(T1w_object) == 0:
			raise ValueError(f"No T1w file found for sub-{subID:02d}, ses-{anatID:02d}.")
		if len(conf_object) == 0:
			raise ValueError(f"No confounds file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")
		if len(reg_object) == 0:
			raise ValueError(f"No TAPAS regressors file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")
		if len(movpar_object) == 0:
			raise ValueError(f"No movement parameters file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")
		if len(out_object) == 0:
			raise ValueError(f"No outliers file found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")
		if len(T1w_to_MNI_object) == 0:
			raise ValueError(f"No T1w to MNI transform file found sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}")
		if len(orig_to_boldref_object) == 0 or len(boldref_to_T1w_object) == 0:
			raise ValueError(f"No boldref to T1w transform files found for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, space-{space}, acq-{acqID}")

		# Check for multiple files (ambiguity)
		if len(log_object) > 1:
			raise ValueError(
				f"Found more than one file:\n{log_object}."
				"\nCheck your LOG folder and logfile import steps."
				)
		
		if len(bold_object) > 1 or len(mask_object) > 1:
			raise ValueError(
				"Found more than one file in one of the following\n:"
				f"{bold_object}, or\n {mask_object}."
				"\n\nPlease check your MRI folder and MRI data import steps."
				)

		if len(conf_object) > 1 or len(movpar_object) > 1 or len(out_object) > 1:
			raise ValueError(
				"Found more than one file in one of the following\n:"
				f"{conf_object},\n{movpar_object}, or\n{out_object}."
				"Check your PHYSIO folder and physiological data import steps."
				)
		
		if len(reg_object) > 1:
			raise ValueError(
				"Found more than one file in one of the following\n:"
				f"{reg_object}."
				"Check your TAPAS preprocessing and regressors estimation steps."
				)
		
		# Warnings for T1 and transform files
		if len(T1w_object) > 1:
			warnings.warn(
				"Multiple anatomical files found. Taking the first one."
			)
		if len(orig_to_boldref_object) > 1 or len(boldref_to_T1w_object) > 1:
			warnings.warn(
                "Multiple transformation files found. Assuming alphabetical order: "
                "orig_to_boldref is the second file, boldref_to_T1w is the first. "
                f"Found: orig={orig_to_boldref_object}, boldref={boldref_to_T1w_object}"
            )
		
		# -------------- 05 Grabing filepaths and Updating --------------
		log_path  = log_object[0].path
		bold_path = bold_object[0].path
		mask_path = mask_object[0].path
		conf_path = conf_object[0].path # selected confonuds
		reg_path  = reg_object[0].path  # only BIOPAC
		movpar_path = movpar_object[0].path # only the trans & rot parameters
		out_path    = out_object[0].path    # motion outliers as detected by fMRIPrep
		T1w_path    = T1w_object[0].path
		T1w_to_MNI_path      = T1w_to_MNI_object[1].path
		orig_to_boldref_path = orig_to_boldref_object[1].path
		boldref_to_T1w_path  = boldref_to_T1w_object[0].path

		# Extract repetition time with PyBIDS methods [sec]
		TR = bold_object[0].get_metadata()['RepetitionTime']
	
	return log_path, bold_path, mask_path, conf_path, reg_path, movpar_path, out_path, T1w_path, T1w_to_MNI_path, orig_to_boldref_path, boldref_to_T1w_path, TR