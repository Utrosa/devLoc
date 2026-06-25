#! /usr/bin/env python
# Time-stamp: <18-05-2026 m.utrosa@bcbl.eu>

# source activate nipypee # loads python libraries
# source ~/.bash_profile # sets SUBJECTS_DIR for freesurfer
# spm is added to matlab path via startup.m

# Import python packages
from pathlib import Path
from nipype import Node, Workflow, IdentityInterface, Function
import nipype.algorithms.modelgen as model
from nipype.interfaces import freesurfer, spm, ants
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Merge

# Import custom-made functions (scripts)
import grabber
from objects_v02 import grab_objects
from designs_v02 import localizer

# -------------------------------------------------------------------------------------------------
# 00. Experiment Parameters
# -------------------------------------------------------------------------------------------------
sub_list = [5]
ses_list = [2, 3, 4, 5, 6, 7]
anat_list = [ses_list[0]]

# Pilot 04 acquisition labels
acqID_list = ["FUNCLOC"]

# Set up project root and define needed folders
homePath   = Path('/home/mutrosa/mutrosa/Documents/devLoc')
tmp_dir    = homePath / "scripts" / "analysis" / "tmp"
out_dir    = homePath / "results"
hrf_dervs  = [0, 0] # using the canonical hrf (without derivatives)
volterra   = False
smoothing  = 2   # Set the Gaussian filter width in mm, default is None
contrasts  = [('localizer', 'T', ['sound', 'silence'], [1, -1])]
MNI        = homePath / "templates" / "tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz" # as in fMRIprep
T1w 	   = homePath / "data_MRI" / "derivatives" / "NORDIC-False" / "sub-05" / "ses-02" / "anat" / "sub-05_ses-02_desc-preproc_T1w.nii.gz"

# -------------------------------------------------------------------------------------------------
# 01. Specify Nodes
# -------------------------------------------------------------------------------------------------
# Infosource: set up a function-free node to iterate over the list of acquisition names.
# The Identity Interface allows to create Nodes that only work with strings (parameters)!
infosource = Node(IdentityInterface(fields = ['subID', 'sesID', 'anatID', 'acqID']),
				  name = "infosource")
infosource.iterables = [('subID', sub_list),
			('sesID', ses_list),
			('anatID', anat_list),
			('acqID', acqID_list)]

# T1w Datasink: create output folder for important outputs in T1w space
datasink_T1w = Node(DataSink(base_directory = str(tmp_dir),
                             container = str(out_dir)),
                name = "datasink_T1w")

# MNI Datasink: create output folder for important outputs in MNI space
datasink_MNI = Node(DataSink(base_directory = str(tmp_dir),
                         	 container = str(out_dir)),
                name = "datasink_MNI")

# Output substitutions: correct all Datasink output folder structures
substitutions = []
subjFolders = [('_acqID_%s_sesID_%s_subID_%s' % (acq, ses, sub),
				'sub-0%s/ses-0%s/acq-%s' % (sub, ses, acq))
               for acq in acqID_list
               for ses in ses_list
               for sub in sub_list]
substitutions.extend(subjFolders)
datasink_T1w.inputs.substitutions = substitutions
datasink_T1w.inputs.substitutions += [('spmT_', 'spmT_space-T1w_'),]
datasink_T1w.inputs.substitutions += [('SPM',   'SPM_space-T1w'),]
datasink_T1w.inputs.substitutions += [('con_',  'con_space-T1w_'),]

datasink_MNI.inputs.substitutions = substitutions
datasink_MNI.inputs.substitutions += [('spmT_', 'spmT_space-MNI_'),]

# Define a Node that extracts filepaths for all files required for the analysis.
infohandle = Node(Function(input_names  = ["subID", "sesID", "anatID", "acqID", "homePath"],
						   output_names = [
						   "log_path", "bold_path", "mask_path", "conf_path",
						   "out_path", "T1w_path", "T1w_toMNI_path", "orig_to_boldref_path",
						   "boldref_to_T1w_path", "TR"
						   				  ], 
						   function = grab_objects),
				name = "infohandle")
infohandle.inputs.homePath = str(homePath)

# Extract information needed to specify a model in the form of a Bunch object.
# Provide info for the Bunch through parsing the event files.
design_bunch = Node(Function(input_names  = ["logfilepath"],
							 output_names = ["design_info"],
							 function = localizer),
					name = "design_bunch")

# Unzip funcional image (preprocessed BOLD).
# This may fail if $SUBJECTS_DIR or $FREESURFER_HOME variables are not set in your profile
unzip = Node(freesurfer.MRIConvert(out_type = 'nii'), name = 'unzip')

# Smoothing
if smoothing is not None:
	smoother = Node(spm.Smooth(fwhm = [smoothing, smoothing, smoothing]),
					name="smooth")

# SpecifyModel: generate SPM-specific godel.
modeler = Node(model.SpecifySPMModel(concatenate_runs = False,
									 input_units  = 'secs',
									 output_units = 'secs',
									 high_pass_filter_cutoff = 128),
			   name = 'modeler')

# Level1Design: generate an SPM design matrix.
designer = Node(spm.Level1Design(bases = {'hrf': {'derivs': hrf_dervs}},
								 timing_units = 'secs',
								 volterra_expansion_order = (2 if volterra else 1)),
				name = 'designer')

# Estimate Model: estimate the parameters of the model.
estimator = Node(spm.EstimateModel(estimation_method = {'Classical': 1}),
				 name = 'estimator')

# Contrast Estimation
contrastor = Node(spm.EstimateContrast(contrasts = contrasts),
				  name = 'contrastor')

# Move data from T1 to MNI space with ANTS. Not necessary if input already in MNI !
# https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.ants.html
warper_toMNI = Node(
		ants.ApplyTransforms(
			dimension = 3,
			interpolation = 'Linear', # Default: Linear
			reference_image = str(MNI),
			invert_transform_flags = [False],  # transform flag set to 0 and as many as transforms (see no. of "transformer" inputs!)
			args = '--float'),
		name = 'warper_toMNI'
	)

warper_toT1w = Node(
		ants.ApplyTransforms(
			dimension = 3,
			interpolation = 'Linear', # Default: Linear
			reference_image = str(T1w),
			invert_transform_flags = [False],  # transform flag set to 0 and as many as transforms (see no. of "transformer" inputs!)
			args = '--float'),
		name = 'warper_toT1w'
	)

# Convert to .nii.gz
zipper_MNI = Node(freesurfer.MRIConvert(out_type = 'niigz'), name = 'zipper_MNI')
zipper_T1w = Node(freesurfer.MRIConvert(out_type = 'niigz'), name = 'zipper_T1w')
# -------------------------------------------------------------------------------------------------
# 02. Connect the Nodes
# -------------------------------------------------------------------------------------------------
l1_localizer = Workflow(name = "l1_localizer")
l1_localizer.base_dir = str(tmp_dir)
l1_localizer.connect([(infosource, infohandle, [("subID", "subID"),
						("sesID", "sesID"),
						("anatID", "anatID"),
						("acqID", "acqID")])])
l1_localizer.connect([(infohandle, design_bunch, [("log_path", "logfilepath")])])
l1_localizer.connect([(infohandle, unzip, [("bold_path", "in_file")])])

if smoothing is not None:
	l1_localizer.connect([(unzip, smoother, [("out_file", "in_files")])])
	l1_localizer.connect([
					(smoother, modeler, [("smoothed_files", "functional_runs")]),
					(infohandle, modeler, [("out_path", "outlier_files"),
										   ("conf_path", "realignment_parameters")]),
					(design_bunch, modeler, [("design_info", "subject_info")]),
					(infohandle, modeler, [("TR", "time_repetition")])
					])
else:
	l1_localizer.connect([
					(unzip, modeler, [("out_file", "functional_runs")]),
					(infohandle, modeler, [
							("out_path", "outlier_files"),
							("conf_path", "realignment_parameters")]),
					(design_bunch, modeler, [("design_info", "subject_info")]),
					(infohandle, modeler, [("TR", "time_repetition")])
					])

l1_localizer.connect([
				(modeler, designer, [("session_info", "session_info")]),
				(infohandle, designer, [("TR", "interscan_interval")])
				])
l1_localizer.connect([
				(designer, estimator, [("spm_mat_file", "spm_mat_file")])
				])
l1_localizer.connect([
				(estimator, contrastor, [("spm_mat_file", "spm_mat_file")]),
				(estimator, contrastor, [("beta_images", "beta_images")]),
				(estimator, contrastor, [("residual_image", "residual_image")]),
				])
l1_localizer.connect([
				(infohandle, warper_toMNI, [('T1w_toMNI_path', 'transforms')]),
				(contrastor, warper_toMNI,  [('spmT_images', 'input_image')])
				])
l1_localizer.connect([
				(infohandle, warper_toT1w, [('boldref_to_T1w_path', 'transforms')]),
				(contrastor, warper_toT1w,  [('spmT_images', 'input_image')])
				])
l1_localizer.connect([
				(contrastor, datasink_T1w, [
					('spm_mat_file', '1stLevel.@spm_mat'),
                               		('spmT_images', '1stLevel.@T'),
                               		('con_images', '1stLevel.@con')])
				])
l1_localizer.connect([
				(warper_toMNI, zipper_MNI, [("output_image", "in_file")])
				])
l1_localizer.connect([
				(warper_toT1w, zipper_T1w, [("output_image", "in_file")])
				])
l1_localizer.connect([
				(zipper_MNI, datasink_MNI, [('out_file', '1stLevel.@T_warped')])
				])
l1_localizer.connect([
				(zipper_T1w, datasink_T1w, [('out_file', '1stLevel.@T_warped')])
				])

# -------------------------------------------------------------------------------------------------
# 03. Visualize the Workflow
# -------------------------------------------------------------------------------------------------
l1_localizer.write_graph(graph2use = 'colored', format = 'png', simple_form = True)

# -------------------------------------------------------------------------------------------------
# 04. Run the Workflow
# -------------------------------------------------------------------------------------------------
res = l1_localizer.run('MultiProc')
