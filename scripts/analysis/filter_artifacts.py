#! /usr/bin/env python
# Time-stamp: <18-06-2026 m.utrosa@bcbl.eu>
"""
Select target artifacts (confounds and motion outliers) for inclusion as 
regressors (realignment parameters) in the first level data analysis.
The artifacts are sourced from fMRIPrep derivatives (timeseries csv) and
TAPAS Physio output.
Note: fMRIPrep estimates confounds from a motion-corrected BOLD images, 
brain mask, mcflirt movement parameters, and a segmentation (source software FSL).

1. CONFOUNDS
fMRIPrep calculates the following confounds:
	- mean global signal,
	- mean tissue class signal,
	- tCompCor, 
	- aCompCor,
	- frame-wise displacement,
	- 6 motion parameters (rot & trans in x, y, z directions),
	- DVARS, and
	- spike regressors.

Defaults confounds are:
	- rot & trans: rigid-body transform parameters indicate how much and how fast you move.
	- csf & wm: the average signal inside cerebrospinal fluid and white-matter masks across time.

Optionally, adds physiological regressors, precomputed with TAPAS
- respiration
- heartbeat
- 02
- CO2

2. MOTION OUTLIERS
Adding only columns for volumes with motion outliers.
Outliers are identified with "1" in the "motion_outlier" columns.
Find the 1 in these columns and extract the volume number (the row number).
"""

# Import python packages
import shutil
import argparse, bids
import pandas as pd
from pathlib import Path
from datetime import datetime

# Import custom-made functions (scripts)
from scripts import grabber

def filter_artifacts(homePath, mriPath, physioPath, subID, sesID, task, denoising, acqIDs, confound_keys, include_biopac=True):
	'''
	Filter artifacts.

	Parameters:
		homePath: Base directory of the project.
		mriPath: The path to the MRI data.
		physioPath: The path to preprocessed physiological regressors files (TAPAS/BIOPAC).
	    subID: Subject identifier.
	    sesID: Session identifier.
		task: Name of the experimental task that the subject was doing. 
		denoising: If True, saving artifacts as denoised with NORDIC.
	    acqIDs: List of acquisition labels in order of data collection.
	    confound_keys: List of confound keys/labels. Rot & Trans must be in FSL order!
	    include_biopac: If True, includes physiological regressors from TAPAS/BIOPAC.

	Side effects:
	    Creates files on disk, including:
	    - txt with confound values (no header)
	    - txt with motion outliers values (no header)
		- txt with motion params (no header)
	    - md summary file including column names of selected artifacts

    '''

	# 00. ---------- Preliminaries ---------- 
	# Define input paths
	homePath = Path(homePath)
	MRILayout = bids.layout.BIDSLayout(mriPath, validate=False, derivatives=True)
	
	# Just a check
	if denoising not in str(mriPath):
		raise ValueError("Are you working on denoised data or not?")

	# Define output paths
	outputPath = homePath / "data_physio" / "artifacts" / f"NORDIC-{denoising}"
	outputPath.mkdir(exist_ok=True, parents=True)

	# Print update to terminal
	print(f"\nReading MRI derivatives data from folder:\n{mriPath}")
	print(f"\nFiltered artifacts will be saved to folder:\n{outputPath}")

	# Specify movement parameters: follow FSL convention for their order
	rot_trans_keys = [
		'rot_x',
		'rot_y',
		'rot_z',
		'trans_x',
		'trans_y',
		'trans_z'
	]
	
	# Take default confounds if none are specified as preferred
	if confound_keys[0] == 'None':
		csf_wm_keys = [
			"csf",
			"csf_derivative1",
			"csf_derivative1_power2", 
			"csf_power2", "white_matter",
			"white_matter_derivative1",
			"white_matter_derivative1_power2",
			"white_matter_power2",
			"csf_wm"
		]
		confounds = rot_trans_keys + csf_wm_keys
	
	# Take specified preferences for confounds
	else:
		confounds = []
		for ck in confound_keys:
			confounds.append(ck)

	# Print update to terminal about the confounds used
	print(f"\nSelecting the following confounds:\n{confounds}") 

	# 01. ---------- Filtering Confounds and Outliers  ---------- 
	outlier_columns = {} # for summary readme file
	for acqID in acqIDs:

		# Grab the correct files
		confound_conf = grabber.define_grabconf(subID, sesID, "timeseries", "tsv", task = task, acquisition = acqID)
		timeseries = grabber.grab_BIDS_object(mriPath, MRILayout, confound_conf)

		if timeseries:
			print(f"\nTimeseries file: {timeseries[0].path}")

			# Load the data
			df_timeseries = pd.read_csv(timeseries[0].path, sep='\t')
			
			# ~~~~~~~~~ CONFOUNDS ~~~~~~~~~
			# Filter the data to select only columns with confound_keys and
			# movement parameters.
			df_movpar = df_timeseries[rot_trans_keys].reset_index(drop=True)
			df_selected = df_timeseries[confounds].reset_index(drop=True)
			
			# Optionally, add BIOPAC regressors as confounds
			if include_biopac:

				# Find physiological confounds from TAPAS
				physioLayout = bids.layout.BIDSLayout(physioPath, validate=False)

				# Grab and load them
				regressor_conf = grabber.define_grabconf(subID, sesID, "regressors", "tsv", acquisition = acqID)
				regressors = grabber.grab_BIDS_object(physioPath, physioLayout, regressor_conf)
				if not regressors:
					raise FileNotFoundError(
						f"No TAPAS/BIOPAC regressor files found in directory: {str(physioPath)} for "
						f"sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acquisitions: {acqIDs}."
					)
				print(f"\nRegressors file: {regressors[0].path}")
				df_regressors = pd.read_csv(regressors[0].path, sep='\t', header=None).reset_index(drop=True)
				
				# Add to dataframe
				print("\nAdding physiological confounds from TAPAS/BIOPAC.\n")
				df_final = df_selected.join(df_regressors)
				df_final = df_final.fillna(0)

			else:
				print("\nNot adding physiological confounds from TAPAS/BIOPAC.\n")
				df_final = df_selected
				df_final = df_final.fillna(0)
			
			# Save movement parameters as a text file.
			movpar_filename = f"sub-{subID:02d}_ses-{sesID:02d}_acq-{acqID}_movpar.txt"
			movpar_out      = outputPath / movpar_filename
			df_movpar.to_csv(movpar_out, sep=' ', header=False, index=False,)
			
			# Save confounds as a text file.		 
			confound_filename = f"sub-{subID:02d}_ses-{sesID:02d}_acq-{acqID}_confounds.txt"
			confounds_out     = outputPath / confound_filename
			df_final.to_csv(confounds_out, sep=' ', header=False, index=False,)

			# ~~~~~~~~~ OUTLIERS ~~~~~~~~~
			# Identify the rows with motion outliers.
			# Note, indexing starts with 0 in pandas and 1 in LibreOffice ;)
			outlier_vols = []
			outlier_cols = []
			for key in df_timeseries.keys():
				if key.startswith('motion_outlier'):
					outlier_cols.append(key)
					for row, value in enumerate(df_timeseries[key]):
						if value == 1:
							outlier_vols.append(row)
			
			# Save column names per acquistion for tracking filtering process
			outlier_columns[acqID] = outlier_cols
			
			# Save outliers as a text file
			outlier_filename = f"sub-{subID:02d}_ses-{sesID:02d}_acq-{acqID}_outliers.txt"
			outlier_out      = outputPath / outlier_filename
			with open(outlier_out, 'w') as f:
				f.write('\n'.join(map(str, outlier_vols)))

			# Print progress report
			print(f"~~~~~~~~ Filtered confounds & artifacts for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID} ~~~~~~~~")
		else:
			raise FileNotFoundError(
				f"Timeseries ({timeseries}).tsv file was not found in {mriPath} for sub-{subID:02d}, ses-{sesID:02d}, task-{task}, acq-{acqID}.")

	# 02. ---------- Save a summary report on which counfounds and outliers were selected.
	readme_file = outputPath / f"sub-{subID:02d}_ses-{sesID:02d}_summary.md"
	with open(readme_file, 'w') as f:
		f.write("# Filer Artifacts - Summary\n")
		f.write(f"- Time and date: {datetime.now()}\n")
		f.write(f"- Acquisitions: {acqIDs}\n")
		f.write(f"- Confounds: {confounds} + BIOPAC {include_biopac}\n")
		f.write(f"- Motion Outliers: {outlier_columns}\n")

	# 03. ---------- Move BIOPAC regressors file to artifacts folder ---------- 
	## Grab the file
	reg_conf   = grabber.define_grabconf(subID, sesID, "regressors", "tsv", acquisition=acqID)
	reg_object = grabber.grab_BIDS_object(physioPath, physioLayout, reg_conf) 
	reg_path   = reg_object[0].path
	regPath    = Path(reg_path)

	## Move it
	dest_path = outputPath / regPath.name
	shutil.copy2(str(regPath), str(dest_path))
	print(f"Moved: {reg_path} -> {dest_path}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument("homePath", type=str)
	parser.add_argument("mriPath", type=str)
	parser.add_argument("physioPath", type=str)
	parser.add_argument("subID", type=int)
	parser.add_argument("sesID", type=int)
	parser.add_argument("task", type=str)
	parser.add_argument("denoising", type=str)
	parser.add_argument("--acqIDs", nargs="+", required=True)
	parser.add_argument("--confound_keys", nargs="+", required=True)
	parser.add_argument(
		"--include_biopac",
		action="store_true",
		help="Includes physiological regressors obtained with TAPAS from BIOPAC data."
		)

	args = parser.parse_args()

	filter_artifacts(
		args.homePath,
		args.mriPath,
		args.physioPath,
		args.subID,
		args.sesID,
		args.task,
		args.denoising,
		args.acqIDs,
		args.confound_keys,
		include_biopac=args.include_biopac
		)