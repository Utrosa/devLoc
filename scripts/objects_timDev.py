#! /usr/bin/env python
# Time-stamp: <2026-05-26 m.utrosa@bcbl.eu>
'''
Grabs objects needed for 1st level GLM analysis in Nipype.
Grabs functional and anatomical files in native space (T1w) !
'''

def grab_objects(subID, sessions, anatID, runs, homePath, mriPath):
	import bids
	import grabber
	from pathlib import Path

	log_paths = []
	bold_paths = []
	mask_paths = []
	conf_paths = []
	out_paths = []
	T1w_paths = []
	T1w_to_MNI_paths = []
	orig_to_boldref_paths = []
	boldref_to_T1w_paths = []
	TRs = []

	for sesID in sessions:
		for runID in runs:
			homePath = Path(homePath)
			mriPath  = Path(mriPath)

			# -------------- 01 Set up layouts -------------- 
			mriLayout = bids.layout.BIDSLayout(mriPath, validate=False, derivatives=True)
			
			logpath   = homePath / "data_logs" / "bids" # CHECK: is this obsolete?
			logLayout = bids.layout.BIDSLayout(logpath, validate=False)

			artpath   = homePath / "data_physio"
			artLayout = bids.layout.BIDSLayout(artpath, validate=False)
			
			# -------------- 02 Configuration -------------- 
			## Log files
			log_conf = grabber.define_grabconf(subID, sesID, "events", "tsv", acquisition = runID)
			
			## Funcional files
			bold_conf   = grabber.define_grabconf(subID, sesID, "bold", "nii.gz", acquisition = runID, space = "T1w") # space = "MNI152NLin2009cAsym"
			bold_object = grabber.grab_BIDS_object(mriPath, mriLayout, bold_conf)
			
			## Outliers and confounds
			mask_conf = grabber.define_grabconf(subID, sesID, "mask", "nii.gz", acquisition = runID, space = "T1w") # space = "MNI152NLin2009cAsym"
			conf_conf = grabber.define_grabconf(subID, sesID, "confounds", "txt", acquisition = runID)
			out_conf  = grabber.define_grabconf(subID, sesID, "outliers",  "txt", acquisition = runID)
			
			## Transform files
			boldref_to_T1w_conf  = grabber.define_grabconf(subID, sesID, "xfm",  "txt")

			## Anatomical files
			T1w_conf        = grabber.define_grabconf(subID, anatID, "T1w",  "nii.gz") # space = "MNI152NLin2009cAsym"
			T1w_to_MNI_conf = grabber.define_grabconf(subID, anatID, "xfm",  "h5")

			# -------------- 03 Grabbing files -------------- 
			log_path = grabber.grab_BIDS_object(logpath, logLayout, log_conf)[0].path
			bold_path  = bold_object[0].path
			mask_path  = grabber.grab_BIDS_object(mriPath, mriLayout, mask_conf)[0].path
			conf_path  = grabber.grab_BIDS_object(artpath, artLayout, conf_conf)[0].path
			out_path   = grabber.grab_BIDS_object(artpath, artLayout, out_conf)[0].path
			T1w_path   = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_conf)[0].path
			T1w_to_MNI_path      = grabber.grab_BIDS_object(mriPath, mriLayout, T1w_to_MNI_conf)[1].path
			orig_to_boldref_path = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)[1].path
			boldref_to_T1w_path  = grabber.grab_BIDS_object(mriPath, mriLayout, boldref_to_T1w_conf)[0].path

			# Extract repetition time with PyBIDS methods [sec]
			TR = bold_object[0].get_metadata()['RepetitionTime']
			
			#  -------------- 04 Append --------------
			log_paths.append(log_path)
			bold_paths.append(bold_path)
			mask_paths.append(mask_path)
			conf_paths.append(conf_path)
			out_paths.append(out_path)
			if len(T1w_paths) < 1:
				T1w_paths.append(T1w_path)
				T1w_to_MNI_paths.append(T1w_to_MNI_path)
			orig_to_boldref_paths.append(orig_to_boldref_path)
			boldref_to_T1w_paths.append(boldref_to_T1w_path)
			TRs.append(TR)

	return log_paths, bold_paths, mask_paths, conf_paths, out_paths, T1w_paths, T1w_to_MNI_paths, orig_to_boldref_paths, boldref_to_T1w_paths, TRs[0]
