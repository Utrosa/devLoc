# Deviance Location Pilot
This is a pilot project to measure neural and behavioral data in response to deviant stimuli.

We conducted a case pilot study. The participant was asked to count how many tones have a deviant pitch from the rest in rhythmic (regular) 7-tone sequences. Some sequences included tones which were timing deviants, meaning that they occured a bit sooner or later than the expected regular tone onset (the beat).

## Project Goals
1. Validate the selected functional 2D EPI sequence with the behavioral paradigm.
2. Test difficulty of the distractor task.

## Summary of Data Analysis Steps 
Data are stored in 3 different folders, depending on their source (MRI scanner, [Expyriment](https://expyriment.org/) software for the behavioral task, and BIOPAC for physiological data):
   - data_logs
   - data_MRI
   - data_physio

### 01 Curation
1. Get the data from the source and save:
	- .dcm files from MRI scanner in data_MRI/sourcedata/dicoms
	- .acq files from BIOPAC in data_physio/sourcedata
	- .txt files from Expyriment in BIDS data_logs/sourcedata
	--> These folders are untouched by next steps to ensure replicable pipeline.
2. Check that data is complete and correctly named. Naming conventions are:
   - .dcm data: sub-{subID:02d}_ses-{sesID:02d}_{project} folders
   - .acq data: sub-{subID:02d}_ses-{sesID:02d}_task-{project}_physio.acq files (localizer vs devLoc)
   - .txt data: include only files from the "bids_output" folder which have ".tsv" extension. All files should be in "sourcedata" folder without subfolders. Task names: localizer, freqDev, timDev.
3. Exclude incomplete data (e.g.: functional scans that were interrupted) and remove any duplicate data. Refer to the laboratory log to guide decisions. The log contains info on execution of the MRI protocol during data acquisition such as errors or modifications.
4. Check onsets in log files.
   - All onsets are correct except for session 02 for frequency counting task (freqDev).
5. Run `checksum.py` to check no files are corrupted.
   - No corrupted files detected in no session.
6. Run `00_pre_import.py` to create sidecar files, needed for the config file.
7. Set up the config file for BIDSifying MRI data with `dcm2bids`. The configuration doesn't have to include the headscouts and the phoenix ZIP report. Validate the config.json file: https://jsonlint.com/.
   - PhaseEncodingDirection is i- (Left-Right; ifmap) & i (Right-Left; fmap)!

### 02 Importing
1. Run: `bash 01a_import_curate.sh`

   **Outputs**:
      - raw MRI data in BIDS
      - background-corrected T1 (mp2rage) images
      - denoised functional and sbref images ([NORDIC](https://github.com/SteenMoeller/NORDIC_Raw/tree/main))
      - removed noise scans from bold and FH sbref scans
      - correctly named and formatted logfiles in (sub-XX/ses-XX/func/)
      - preprocessed physio data and physiological noise regressors (TAPAS) per session, subject, and functional sequence. Note, there's an iteration-over-aquisitions-of-same-task option in tapas.m)

   **Warnings**:
   * WARNING | Chris Rorden's dcm2niiX version v1.0.20250505  GCC10.2.1 x86-64 (64-bit Linux)
   * Warning: 4D Siemens XA images should be exported as enhanced not classic DICOM. Slice times and other properties may be inaccurate.
   * Warning: X does not support locale en_US.UTF-8 (while running nordic in MATLAB)
2. Exclude incomplete data (e.g.: functional scans that were interrupted) and remove any duplicate data. Refer to the laboratory log to guide decisions. The log contains info on execution of the MRI protocol during data acquisition such as errors or modifications.
   - Shorten functional and phasic data of task-localizer_acq-FUNCLOC scan, collected in session 03. Original nvols was 511. Shortened to 395 (accounting for noise scan).
   - Split physio data for session 02 and 06 into "timDev" and "localizer" files
3. Visually inspect images by running `bash 01b_visualize.sh`.
4. Add task stimuli, `dataset_description`, and `README` files to data_MRI/sourcedata/raw/.
5. Run [BIDS Validator](http://bids.neuroimaging.io/tools/validator.html) on the dataset to ensure compliance to [the latest BIDS specification](https://bids-specification.readthedocs.io/en/stable/).

### 03 Preprocessing
1. Run: `bash 02_fMRIprep.sh`
   Outputs:
   - preprocessed fMRI data
2. Delete temporary cache and work directories once preprocessing is successful.

### 04 1st Level Analysis: GLM
1. Download the subcortical atlas and MNI template.
   - [Sitek's in-vivo subcortical atlas](https://github.com/sitek/subcortical-auditory-atlas/tree/master/atlases)
   - [MNI template from Template Flow](https://www.templateflow.org/archive/)
2. Run `python resample_atlas.py` to:
   - resample Sitek's in-vivo atlas to the resolution of the MNI template used in preprocessing and data analyses, and then to T1w native space (for subcortical ROIS)
   - resample Freesurfer's reconall atlas to T1w native and MNI spaces (for cortical ROIS)
3. Run `bash 03a_filer_artifacts.sh`

   Note, raw physiological data, collected with BIOPAC, is independent from NORDIC denoising steps, while artifact physiological data is not. The confounds text file, created by `filter_artifacts.py`, contains the NORDIC-independent physiological confounds (from [RETROICOR model](https://doi.org/10.1002/1522-2594(200007)44:1%3C162::AID-MRM23%3E3.0.CO;2-E)) and selected confounds from fMRIPrep timeseries file (translations and rotations).

   Outputs:
   - the selected confounds per volume (physiological artifacts and selected confounds from fMRIPrep - FSL mcflirt)
   - motion outliers (as caluculated by fMRIPrep - FSL mcflirt)
   - motion parameters (translations & rotations)
4. Run `bash 03b_analyze_task-timDev.sh` or `bash 03b_analyze_task-localizer.sh`

   Outputs per subject, session, task, acquisition, and optionally, run:
      - betas 
      - residuals
      - contrasts
      - SPM design

5. Run `python resample_outputs.py` to:
   Before extracting the timeseries from the ROIs: resample the outputs from the 1st level analysis (betas/contrasts/t-values) from restricted FoV space of functional scans to native space of the T1w image using the "from_boldref_to_T1w" transformation file. This ensures that the outputs and the atlas have the same shape and affines.

6. Descriptive plotting
   Returns:
      - 

#### 05. 2nd Level Analysis
5. Run `.py`

   Performs non-parametric tests per voxel of an ROI or per run (averaging ROI voxels)
   Outputs averaged arrays (across experimental runs and sessions) per specifed region of interests (ROIs) and plots them.

# License
This project is licensed under the terms of the MIT License.