#!/usr/bin/env python
# Time-stamp: <04-05-2026 m.utrosa@bcbl.eu>
"""
Prerequisites:
- For python: phys2bids; bioread
- For matlab: Signal Processing; Statistics and Machine Learning

For splitting raw physiological data into separate files per task.
"""

# Import python packages
import bids, subprocess, shutil, sys, warnings
from pathlib import Path

# Import custom-made functions
from scripts import grabber

def split_PHYSIO(subID, sesID, project, task, homePath):
	"""
	Import raw physiological data and organize them into a BIDS-compliant structure.

	Parameters:
	    subID: Subject identifier.
	    sesID: Session identifier.
	    project: Project name as entered in the Siemens computer.
	    task: Name of the experimental task that the subject was doing.
	    homePath: Base directory of the project.
	
	Raises:
	    FileNotFoundError if required BIOPAC files are missing.
	    ValueError if .acq files are not found (due to incorrect name).
	"""

	# 1. Physio data import -------------------------------------------------------------
	homePath   = Path(homePath)
	physioPath = homePath / "data_physio" / "sourcedata"
	if not physioPath.exists():
		raise FileNotFoundError(f"Physio folder not found: {physioPath}")

	# Find .acq files following BIDS-compliant naming convention
	physioLayout = bids.layout.BIDSLayout(physioPath, validate=False)
	physio_conf  = grabber.define_grabconf(subID, sesID, "physio", "acq", task=task)
	acq_object   = grabber.grab_BIDS_object(physioPath, physioLayout, physio_conf)
	if len(acq_object) == 0:
		files = "\n".join(item.name for item in physioPath.iterdir())
		raise ValueError(f"\nUnexpected physio name format: \n{files}")

	# Check the amount of .acq files collected. Expecting one file per session.
	acq_path = Path(acq_object[0].path)
	acqFile  = acq_path.name
	if len(acq_object) > 1:
		warnings.warn(
			f"\nMultiple BIOPAC files found ({len(acq_object)})"
			f"For subject {subID}, session {sesID}, task {task} using {acq_object[1].name}.")
	print(f"\nPhysio folder: {physioPath}")
	print(f"\nGrabbing: {acqFile}")

	# 2. Create a temporary folder to store TAPAS-compatible data -----------------------
	rawPath = homePath / "data_physio" / "raw" 
	rawPath.mkdir(exist_ok=True, parents=True)
	tmpPath = homePath / "data_physio" / "tmp"
	tmpPath.mkdir(exist_ok=True, parents=True)

	# 3. Transform .acq files to BIDS-compliant tsv -------------------------------------
	# Splitting by number of volumens (main func + FH func scans) and TR (1.51)
	ntp_values = [511, 2, 511, 2, 511, 2, 511, 2, 396, 2]
	tr_values = [1.51, 1.51, 1.51, 1.51, 1.51, 1.51, 1.51, 1.51, 1.51, 1.51]
	bids_cmd = [
			"phys2bids",
			"-in", acqFile,
			"-indir", str(physioPath),
			"-chtrig", str(1),
			"-ntp",  *map(str, ntp_values),
			"-tr", *map(str, tr_values),
			"-thr", str(1),
			"-outdir", str(tmpPath),
	]

	print(f"\nRunning: {' '.join(bids_cmd)}")
	subprocess.run(bids_cmd, check=True)
	
if __name__ == "__main__":
	subID, sesID, project, task, homePath = (
		int(sys.argv[1]),
		int(sys.argv[2]),
		sys.argv[3],
		sys.argv[4],
		sys.argv[5]
		)
	split_PHYSIO(subID, sesID, project, task, homePath)