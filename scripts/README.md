# Deviance Location Pilot
### Fixed Effects Analyses
1. Prepare the atlases.
   - Download the [Sitek's in-vivo subcortical atlas](https://github.com/sitek/subcortical-auditory-atlas/tree/master/atlases)
   - Use freesurfer's utility function `[mri_convert](https://surfer.nmr.mgh.harvard.edu/fswiki/mri_convert)` to covert the recon-all parcellation atlases from MGH to NifTi format. Note, recon-all will parcellate the individual subject’s brain according to two atlases: the Desikan-Killiany atlas (aparc+aseg.mgz) and the Destrieux atlas (aparc.a2009s+aseg.mgz). The Destrieux atlas contains more parcellation.

2. Resample the atlases to the shape of the results (e.g.: SPM t-values or beta images). Run the appropriate section of `python resample_atlas.py`. Please do not run the former script as a whole!

3. Run `bash 03a_filer_artifacts.sh`
	Note, raw physiological data, collected with BIOPAC, is independent from NORDIC denoising steps, while artifact physiological data is not. The confounds text file, created by `filter_artifacts.py`, contains the NORDIC-independent physiological confounds (from [RETROICOR model](https://doi.org/10.1002/1522-2594(200007)44:1%3C162::AID-MRM23%3E3.0.CO;2-E)) and selected confounds from fMRIPrep timeseries file (translations and rotations).

   Outputs:
	- the selected confounds per volume (physiological artifacts and selected confounds from fMRIPrep - FSL mcflirt)
	- motion outliers (as caluculated by fMRIPrep - FSL mcflirt)
	- motion parameters (translations & rotations)

4. Run `python analysis_GLM.py`
   Outputs per subject, session, task, aquisition, and optionally, run:
	- beta images from the SPM design matrix
	- residual images
	- SPM design image

5. Run `average_betas.py`
Before extracting the timeseries from the ROIs: resample the atlas and resulting nifti files (beta images or t-values) to T1 space using the "from_boldref_to_T1w" transformation file. This ensures that the resulting image and the atlas have the same shape and affines.

   Outputs averaged arrays (across experimental blocks and sessions) per specifed region of interests (ROIs) and plots them.

