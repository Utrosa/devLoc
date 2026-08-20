#!/usr/bin/env python
import bioread, glob, os, shutil

# Ensure that you have all the BIDS-compliant folders. Run this command only once.
os.system("dcm2bids_scaffold")

# 1. Adjust parameters for your subject/session
# Activate the correct conda environment: source activate localizer_fMRI.
subID    = 5
sesID    = 7
project  = "devLoc"
homePath = '/home/mutrosa/mutrosa/Documents/devLoc/data_MRI/'

# 2. Generate sidecar files
# These are needed for setting up the conf.json and checking the completeness of data.
# DICOM directory(ies) or archive(s)
helper_command =  f"dcm2bids_helper -d {homePath}sourcedata/dicoms/sub-{subID:02d}_ses-{sesID:02d}_{project} "

# Output directory
helper_command += f"-o {homePath}sidecars/sub-{subID:02d}/ses-{sesID:02d}"
os.system(helper_command)
