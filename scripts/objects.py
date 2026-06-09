#! /usr/bin/env python
# Time-stamp: <2026-03-06 m.utrosa@bcbl.eu>
'''
Grabs objects needed for 1st level GLM analysis in Nipype.
Grabs functional and anatomical files in the specified space.

If sessions is a list of sessions, returns a list of objects for the
specified sessions. Same for runs.

TODO: add option in which the code doesn't fail if a certain acq was not collected in a session ... make more modular.
TODO: test that it works as expected with runs and MNI
'''

def grab_objects(subID, sesID, anatID, homePath, mriPath, space, acqID, runs):
	"""
	Locate functional and anatomical objects and return a flat list of filepaths of those objects
	based on the specified subject, sessions, acquisitions, and runs.

	Args:
		subID (int): The subject identifier.
		sesID (str): Session identifier.
		anatID (int): The identifier of the session in which the anatomical image was obtained.
		homePath (str): The base directory path. Defaults to None.
		mriPath (str): The path to the MRI data. Defaults to None.
		space (str): The target space for the files (e.g., MNI152NLin2009cAsym or T1w).
		acqID (str): Acquisition name.
		runs (bool, str): If False, no run identifier is used to locate files.

	Returns:
		lists: Flat lists of paths to objects corresponding to the specified parameters.

	Note:
		Ensure that sessions, acquisitions, and runs match.
	"""
	import bids
	import grabber
	from pathlib import Path

	if runs:					
		homePath = Path(homePath)
		mriPath  = Path(mriPath)

		# -------------- 01 Set up layouts -------------- 
		logpath   = homePath / "data_logs" / "bids"
		logLayout = bids.layout.BIDSLayout(logpath, validate=False)
		mriLayout = bids.layout.BIDSLayout(mriPath, validate=False, derivatives=True)
		artpath   = homePath / "data_physio"
		artLayout = bids.layout.BIDSLayout(artpath, validate=False)
		
		# -------------- 02 Configuration -------------- 
		log_conf  = grabber.define_grabconf(subID, sesID, "events", "tsv", acquisition = acqID, run = f"{runID:02d}")
		bold_conf = grabber.define_grabconf(subID, sesID, "bold", "nii.gz", acquisition = acqID, run = f"{runID:02d}", space = space)
		mask_conf = grabber.define_grabconf(subID, sesID, "mask", "nii.gz", acquisition = acqID, run = f"{runID:02d}", space = space)
		conf_conf = grabber.define_grabconf(subID, sesID, "confounds", "txt", acquisition = acqID, run = f"{runID:02d}")
		out_conf  = grabber.define_grabconf(subID, sesID, "outliers",  "txt", acquisition = acqID, run = f"{runID:02d}")
		T1w_conf  = grabber.define_grabconf(subID, anatID, "T1w",  "nii.gz")
		T1w_to_MNI_conf = grabber.define_grabconf(subID, anatID, "xfm",  "h5")
		boldref_to_T1w_conf  = grabber.define_grabconf(subID, sesID, "xfm",  "txt")

		# -------------- 03 Grabbing files -------------- 
		log_path    = grabber.grab_BIDS_object(logpath, logLayout, log_conf)[0].path
		bold_object = grabber.grab_BIDS_object(mriPath, mriLayout, bold_conf)
		bold_path   = bold_object[0].path
		mask_path   = grabber.grab_BIDS_object(mriPath, mriLayout, mask_conf)[0].path
		conf_path   = grabber.grab_BIDS_object(artpath, artLayout, conf_conf)[0].path
		out_path    = grabber.grab_BIDS_object(artpath, artLayout, out_conf)[0].path
		T1w_path    = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_conf)[0].path
		T1w_to_MNI_path      = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_to_MNI_conf)[1].path
		orig_to_boldref_path = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)[1].path
		boldref_to_T1w_path  = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)[0].path

		# Extract repetition time with PyBIDS methods [sec]
		TR = bold_object[0].get_metadata()['RepetitionTime']
			
	else:
		homePath = Path(homePath)
		mriPath  = Path(mriPath)

		# -------------- 01 Set up layouts -------------- 
		logpath   = homePath / "data_logs" / "bids"
		logLayout = bids.layout.BIDSLayout(logpath, validate=False)
		mriLayout = bids.layout.BIDSLayout(mriPath, validate=False, derivatives=True)
		artpath   = homePath / "data_physio"
		artLayout = bids.layout.BIDSLayout(artpath, validate=False)
		
		# -------------- 02 Configuration -------------- 
		log_conf  = grabber.define_grabconf(subID, sesID, "events", "tsv", acquisition = acqID)
		bold_conf = grabber.define_grabconf(subID, sesID, "bold", "nii.gz", acquisition = acqID, space = space)
		mask_conf = grabber.define_grabconf(subID, sesID, "mask", "nii.gz", acquisition = acqID, space = space)
		conf_conf = grabber.define_grabconf(subID, sesID, "confounds", "txt", acquisition = acqID)
		out_conf  = grabber.define_grabconf(subID, sesID, "outliers",  "txt", acquisition = acqID)
		T1w_conf  = grabber.define_grabconf(subID, anatID, "T1w",  "nii.gz")
		T1w_to_MNI_conf = grabber.define_grabconf(subID, anatID, "xfm",  "h5")
		boldref_to_T1w_conf  = grabber.define_grabconf(subID, sesID, "xfm",  "txt")

		# -------------- 03 Grabbing files -------------- 
		log_path    = grabber.grab_BIDS_object(logpath, logLayout, log_conf)[0].path
		bold_object = grabber.grab_BIDS_object(mriPath, mriLayout, bold_conf)
		bold_path   = bold_object[0].path
		mask_path   = grabber.grab_BIDS_object(mriPath, mriLayout, mask_conf)[0].path
		conf_path   = grabber.grab_BIDS_object(artpath, artLayout, conf_conf)[0].path
		out_path    = grabber.grab_BIDS_object(artpath, artLayout, out_conf)[0].path
		T1w_path    = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_conf)[0].path
		T1w_to_MNI_path      = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_to_MNI_conf)[1].path
		orig_to_boldref_path = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)[1].path
		boldref_to_T1w_path  = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)[0].path

		# Extract repetition time with PyBIDS methods [sec]
		TR = bold_object[0].get_metadata()['RepetitionTime']
	
	return log_path, bold_path, mask_path, conf_path, out_path, T1w_path, T1w_to_MNI_path, orig_to_boldref_path, boldref_to_T1w_path, TR