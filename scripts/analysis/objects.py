#! /usr/bin/env python
# Time-stamp: <2026-15-06 m.utrosa@bcbl.eu>
'''
Grabs objects needed for other scripts.
Grabs functional and anatomical files in the specified space.
'''
def check_object(obj, name, sub, ses, extra_params="", warning_only=False):
    """
    Raises a ValueError for missing or ambiguous files.

    Args:
        obj (list): The list of file objects to validate.
        name (str): A descriptive name for the file type.
        sub (int): The subject identifier.
        ses (str): The session identifier.
        extra_params (str): Additional context parameters to append to the message.
        warning_only (bool): If True, issue a warning instead of raising an error.
    """
    import warnings
    count = len(obj)

    # Check for missing values
    if count == 0:
        msg = f"No {name} found for sub-{sub:02d}, ses-{ses:02d}, {extra_params}."
        if warning_only:
            warnings.warn(msg)
        else:
            raise ValueError(msg)
	
	# Check which files are found as a group if more than one file found   
    elif count > 1:
        msg = f"Found more than one {name}:\n{obj}.\nPlease verify your file-grabbing inputs."
        warnings.warn(msg)

def grab_objects(subID, sesID, anatID, homePath, mriPath, artPath, space, task, acq=None, run=None):
	"""
	Locate functional and anatomical objects and returns a tuple of filepaths and TR based on the 
	specified subject and session. Optionally, you can specify the acquisition and run.

	Args:
		subID (int): The subject identifier.
		sesID (str): Session identifier.
		anatID (int): The identifier of the session in which the anatomical image was obtained.
		homePath (str): The base directory path.
		mriPath (str): The path to the MRI data.
		artPath (str): ThE path to the preprocessed physiological data.
		space (str): The target space for the files (e.g., MNI152NLin2009cAsym or T1w).
		task (str): Name of the experimental task that the subject was doing.

	Òptional args:
		acq (None/str): If None, no acquisition name is used to locate files. If str (True),
						that acquisition identifier is used to find files.
		run (None/str): If None, no run identifier is used to locate files. If str (True),
						that run identifier is used to find files.

	Returns:
		tuple: strings (filepaths) and TR (float).

	Note:
		Ensure that session, acquisition, and run match an existing scenario.
	"""
	
	import bids
	import grabber
	import warnings
	from pathlib import Path	

	# Initialize paths
	homePath = Path(homePath)
	mriPath  = Path(mriPath)
	artPath  = Path(artPath)

	# -------------- 01 Set up layouts -------------- 
	logpath   = homePath / "data_logs" / "bids"
	logLayout = bids.layout.BIDSLayout(logpath, validate=False)
	mriLayout = bids.layout.BIDSLayout(mriPath, validate=False)
	artLayout = bids.layout.BIDSLayout(artPath, validate=False)

	# Determine if run and functional acquisition identifiers are used
	runID = run	if bool(run) else None
	acqID = acq if bool(acq) else None
	
	# -------------- 02 Configuration -------------- 
	log_conf = grabber.define_grabconf(subID, sesID, "events", "tsv", task=task, acquisition=acqID, run=runID)
	bold_conf = grabber.define_grabconf(subID, sesID, "bold", "nii.gz", task=task, acquisition=acqID, run=runID, space=space)
	mask_conf = grabber.define_grabconf(subID, sesID, "mask", "nii.gz", task=task, acquisition=acqID, run=runID, space=space)
	conf_conf = grabber.define_grabconf(subID, sesID, "confounds", "txt", task=task, acquisition=acqID, run=runID)	
	reg_conf = grabber.define_grabconf(subID, sesID, "regressors", "tsv", task=task, acquisition=acqID, run=runID)
	movpar_conf = grabber.define_grabconf(subID, sesID, "movpar", "txt", task=task, acquisition=acqID, run=runID)
	out_conf = grabber.define_grabconf(subID, sesID, "outliers",  "txt", task=task, acquisition=acqID, run=runID)
	T1w_conf = grabber.define_grabconf(subID, anatID, "T1w",  "nii.gz")
	h5_trans_conf = grabber.define_grabconf(subID, anatID, "xfm",  "h5")
	txt_trans_conf = grabber.define_grabconf(subID, anatID, "xfm",  "txt", task=task, acquisition=acqID)

	# -------------- 03 Grabbing files --------------
	log_object = grabber.grab_BIDS_object(logpath, logLayout, log_conf)
	bold_object = grabber.grab_BIDS_object(mriPath, mriLayout, bold_conf)
	mask_object = grabber.grab_BIDS_object(mriPath, mriLayout, mask_conf)
	conf_object = grabber.grab_BIDS_object(artPath, artLayout, conf_conf) # selected confounds
	reg_object = grabber.grab_BIDS_object(artPath, artLayout, reg_conf)   # only BIOPAC
	movpar_object = grabber.grab_BIDS_object(artPath, artLayout, movpar_conf) # only the trans & rot parameters
	out_object = grabber.grab_BIDS_object(artPath, artLayout, out_conf) # motion outliers as detected by fMRIPrep
	T1w_object = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_conf)
	T1w_to_MNI_object = grabber.grab_BIDS_object(mriPath, mriLayout, h5_trans_conf)
	func_trans_object = grabber.grab_BIDS_object(mriPath, mriLayout, txt_trans_conf)

	# -------------- 04 Verification & Warnings --------------
	space_str = f", space-{space}" if space else ""
	task_str = f", task-{task}" if task else ""
	acq_str = f", acq-{acqID}" if bool(acq) else ""
	run_str = f", run-{runID:02d}" if bool(run) else ""
	extra_str = f"{task_str}, {space_str}, {acq_str}, {run_str}"

	# Missing file checks
	check_object(log_object, "log file", subID, sesID, extra_str, warning_only=True)
	check_object(bold_object, "bold file", subID, sesID, extra_str, warning_only=True)
	check_object(mask_object, "mask file", subID, sesID, extra_str, warning_only=True)
	check_object(T1w_object, "T1w file", subID, anatID)
	check_object(conf_object, "confounds file", subID, sesID, extra_str)
	check_object(reg_object, "TAPAS regressors file", subID, sesID, f"{acq_str}, {run_str}", warning_only=True)
	check_object(movpar_object, "movement parameters file", subID, sesID, extra_str)
	check_object(out_object, "outliers file", subID, sesID, extra_str)
	check_object(T1w_to_MNI_object, "T1w to MNI transform file", subID, sesID, extra_str)

	# Transform files
	if len(func_trans_object) == 0:
		raise ValueError(
			f"No orig_to_boldref or boldref_to_T1w transform files found for sub-{subID:02d}, ses-{sesID:02d}{extra_str}"
		)

	# Warnings for multiple files
	if len(T1w_object) > 1:
		warnings.warn(f"Multiple anatomical files found: {[Path(to).name for to in T1w_object]}")

	if len(func_trans_object) > 1:
		names_orig = [Path(otbo).name for otbo in func_trans_object]
		names_bold = [Path(btto).name for btto in func_trans_object]
		warnings.warn(
			f"Multiple transformation files found: \n"
			f" * orig_to_boldref: {names_orig} \n\n"
			f" * boldref_to_T1w: {names_bold}"
		)

	# -------------- 05 Grabing filepaths and Updating --------------
	log_paths    = [lo.path for lo in log_object]
	bold_paths   = [bo.path for bo in bold_object]
	mask_paths   = [mo.path for mo in mask_object]
	conf_paths   = [co.path for co in conf_object] # selected confounds
	movpar_paths = [mpo.path for mpo in movpar_object]
	out_paths    = [oo.path for oo in out_object]  # motion outliers as detected by fMRIPrep

	# Select the anatomical file and print selection to terminal
	T1w_path = T1w_object[0].path
	warnings.warn(
		f"\nThe selected space for the analysis is: {space}. "
		f"\nThe anatomical file selected is: {Path(T1w_object[0]).name}."
	)

	# Regressors tsv file only exists when including BIOPAC regressors
	if len(reg_object):
		reg_paths = [ro.path for ro in reg_object]
	else:
		reg_paths = []

	# Transformation paths
	T1w_to_MNI_path = []
	for ttMNI in T1w_to_MNI_object:
		if "from-T1w_to-MNI" in str(ttMNI):
			T1w_to_MNI_path = ttMNI.path
			print(f"\nFor from-T1w_to-MNI selected: {Path(ttMNI).name}")

	orig_to_boldref_paths = []
	for otbo in func_trans_object:
		if "from-orig_to-boldref" in str(otbo):
			orig_to_boldref_paths.append(otbo.path)
			print(f"\nFor from-orig_to-boldref selected: {Path(otbo).name}")

	boldref_to_T1w_paths = []
	for btto in func_trans_object:
		if "from-boldref_to-T1w" in str(btto):
			boldref_to_T1w_paths.append(btto.path)
			print(f"\nFor from-boldref_to-T1w selected: {Path(btto).name}")

	# Extract repetition time with PyBIDS methods [sec]
	TRs = [bo.get_metadata()['RepetitionTime'] for bo in bold_object]
		
	return log_paths, bold_paths, mask_paths, conf_paths, reg_paths, movpar_paths, out_paths, T1w_path, T1w_to_MNI_path, orig_to_boldref_paths, boldref_to_T1w_paths, TRs