# Deviance Location Pilot
This is a pilot project to measure neural and behavioral data in response to deviant stimuli.

We conducted a case pilot study. The participant was asked to count how many tones have a deviant pitch from the rest in rhythmic (regular) 7-tone sequences. Some sequences included tones which were timing deviants, meaning that they occured a bit sooner or later than the expected regular tone onset (the beat).

# Usage
[add examples]

## Project Goals
1. Validate the selected functional 2D EPI sequence with the behavioral paradigm.
2. Test difficulty of the distractor task.

## Summary of Data Analysis Steps 
Data are stored in 3 different folders, depending on their source (MRI scanner, [Expyriment](https://expyriment.org/) software for the behavioral task, and BIOPAC for physiological data):
- data_logs
- data_MRI
- data_physio

### 01 Curation
1. Get the data from the source and save :
	.dcm files from MRI scanner in data_MRI/sourcedata/dicoms
	.acq files from BIOPAC in data_physio/sourcedata
	.txt files from Expyriment in BIDS data_logs/sourcedata
	--> These folders are untouched by next steps to ensure replicable pipeline.
2. Check that data is complete and correctly named. Naming conventions are:
   .dcm data: sub-{subID:02d}_ses-{sesID:02d}_{project} folders
   .acq data: sub-{subID:02d}_ses-{sesID:02d}_task-{project}_physio.acq files (localizer vs devLoc)
   .txt data: include only files from the "bids_output" folder which have ".tsv" extension. All files should be in "sourcedata" folder without subfolders. Task names: localizer, freqDev, timDev.
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
   - preprocessed physio data and physiological noise regressors (TAPAS) per session, subject, and functional sequence (note, there's an iteration-over-aquisitions-of-same-task option in tapas.m)

**Warnings**:
* WARNING | Chris Rorden's dcm2niiX version v1.0.20250505  GCC10.2.1 x86-64 (64-bit Linux)
* Warning: 4D Siemens XA images should be exported as enhanced not classic DICOM. Slice times and other properties may be inaccurate.
* Warning: X does not support locale en_US.UTF-8 (while running nordic in MATLAB)

2. Exclude incomplete data (e.g.: functional scans that were interrupted) and remove any duplicate data. Refer to the laboratory log to guide decisions. The log contains info on execution of the MRI protocol during data acquisition such as errors or modifications.
- DONE: Shorten functional and phasic data of task-localizer_acq-FUNCLOC scan, collected in session 03. Original nvols was 511. Shortened to 395 (accounting for noise scan).
- TODO: Split physio data for session 02 and 06 into "timDev" and "localizer" files

3. Visually inspect images by running `bash 01b_visualize.sh`.
4. Add task stimuli, `dataset_description`, and `README` files to data_MRI/sourcedata/raw/.
5. Run [BIDS Validator](http://bids.neuroimaging.io/tools/validator.html) on the dataset to ensure compliance to [the latest BIDS specification](https://bids-specification.readthedocs.io/en/stable/).

### 03 Importing & Preprocessing
1. Run: `bash 02_fMRIprep.sh`
   Outputs:
   - preprocessed fMRI data
2. Delete temporary cache and work directories once preprocessing is successful.

### 04 Functional Localizer Analysis
1. Download the subcortical atlas and MNI template.
   - [Sitek's in-vivo subcortical atlas](https://github.com/sitek/subcortical-auditory-atlas/tree/master/atlases)
   - [MNI template from Template Flow](https://www.templateflow.org/archive/)
2. Resample Sitek's atlas to MNI or T1w resolution by running `python resample_atlas.py`.
3. Run `bash 03a_filer_artifacts.sh`
   Outputs:
	- modified counfounds and outliers (physiological artifacts)
4. Run `bash 03b_analyze_localizer.sh`
   Outputs:

TODO: explain the resampling from boldref to T1 per func scan & results images (betas, contrasts).
TODO: explain the plotting of betas and contrasts with 2ndLevelAnalysis (non-parametric tests per voxel of an ROI and per run (averaging ROI voxels)).
TODO: list other steps ...


# License
This project is licensed under the terms of the MIT License.



