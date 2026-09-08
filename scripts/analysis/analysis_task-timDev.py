#! /usr/bin/env python
# Time-stamp: <07-09-2026 m.utrosa@bcbl.eu>
'''
fMRI: GLM model fitting with fixed effects

Tutorial
https://nipype.readthedocs.io/en/latest/users/examples/fmri_nipy_glm.html
'''
# Import prerequisites from python and Nipype
from pathlib import Path
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.io import DataSink
from nipype import Workflow, Function, IdentityInterface
import nipype.interfaces.spm as spm  # spm
import nipype.pipeline.engine as pe  # pypeline engine
import nipype.algorithms.modelgen as model  # model specification
from nipype.algorithms.rapidart import ArtifactDetect # artifact detection

# Import custom-made functions (scripts)
import grabber
import config as c
from objects import grab_objects
from designs import timfreqDev
from utils import add_nuisance

# -------------------------------------------------------------------------------------------------
# 01. Specify helper nodes
# -------------------------------------------------------------------------------------------------
# Infosource: set up a function-free node to iterate over the list of acquisition names.
# The Identity Interface allows to create Nodes that only work with strings (parameters)!
infosource = pe.Node(
    IdentityInterface(fields = ['subID', 'sesID', 'acqID']),
	name = "infosource"
)
infosource.iterables = [
    ('subID', c.subIDs),
	('sesID', c.sesIDs),
	('acqID', c.acqIDs)
]

# T1w Datasink: create output folder for important outputs in T1w space
datasink_T1w = pe.Node(
    DataSink(
        base_directory = str(c.workDir),
        container = str(c.outDir)
    ),
    name = "datasink_T1w"
)

# Output substitutions: correct all Datasink output folder structures
substitutions = []
subjFolders = [('_acqID_%s_sesID_%s_subID_%s' % (acq, ses, sub),
				'sub-0%s/ses-0%s/acq-%s' % (sub, ses, acq))
               for acq in c.acqIDs
               for ses in c.sesIDs
               for sub in c.subIDs]
substitutions.extend(subjFolders)
datasink_T1w.inputs.substitutions = substitutions
datasink_T1w.inputs.substitutions += [('beta_', 'beta_space-boldref_'),]
datasink_T1w.inputs.substitutions += [('con_',  'con_space-boldref_'),]
datasink_T1w.inputs.substitutions += [('spmT_',  'spmT_space-boldref_'),]

# Define a Node that extracts filepaths for all files required for the analysis
infohandle = pe.Node(
    Function(
        input_names  = [
            "subID",
            "sesID", 
            "anatID", 
            "homePath", 
            "mriPath",
            "artPath",
            "space",
            "acqID", 
            "task",
            "run"
        ],
        output_names = [
            "log_path",
            "bold_path", 
            "mask_path", 
            "conf_path",
            "reg_path", # empty list (no path to a regressors file) if BIOPAC off (no TAPAS)
            "movpar_path",
            "out_path", 
            "T1w_path", 
            "T1w_to_MNI_path",
            "orig_to_boldref_path",
            "boldref_to_T1w_path", 
            "TR"
        ],
        function = grab_objects),
name = "infohandle"
)
infohandle.inputs.anatID   = c.anatID
infohandle.inputs.homePath = str(c.homePath)
infohandle.inputs.mriPath  = str(c.mriPath)
infohandle.inputs.artPath  = str(c.artPath)
infohandle.inputs.space    = c.space
infohandle.inputs.task     = c.task
infohandle.inputs.run      = False # TODO: check this?

# -------------------------------------------------------------------------------------------------
# 02. Additional preprocessing nodes: smoothing and outlier detection
# -------------------------------------------------------------------------------------------------
if c.smoothing:
    smoother = pe.Node(interface=spm.Smooth(), name="smoother")
    smoother.inputs.fwhm = c.smoothing # TODO: Check whether this has to be a list?  

if c.artDetect:
    # Using intensity and motion parameters to infer parameters
    # https://nipype.readthedocs.io/en/latest/api/generated/nipype.algorithms.rapidart.html
    art_detect = pe.Node(
        ArtifactDetect(), # performs artifact detection on functional images
        name = "art_detect"
    )
    art_detect.inputs.parameter_source = "FSL" # fMRIPrep uses FSL MCFLIRT to estimate confounds
    art_detect.inputs.mask_type = "file"
    art_detect.inputs.save_plot = True # Save plots containing outliers

    # art_detect.inputs.norm_threshold = 1 # Default from documentation's example
    # art_detect.inputs.zintensity_threshold = 3 # Default from documentation's example
    art_detect.inputs.rotation_threxxxxshold    = c.rot_thresh
    art_detect.inputs.translation_threshold = c.trans_thresh

    # Deterimne which differences to use for outlier detection: Motion and Intensity parameters
    art_detect.inputs.use_differences = [True, False]

# -------------------------------------------------------------------------------------------------
# 03. Specify 1st-level model parameters
# -------------------------------------------------------------------------------------------------
# Get the information about the experimental paradigm to create an SPM design matrix.
# Construct a list of objects (each object should contain data for all runs of that session)
# Create a Bunch object by parsing all event files: timDev & freqDev are separate Bunch objects.
bunch_log = pe.Node(
    Function(
        input_names = ["time_log", "time_groups", "time_pool", "time_binary"],
        output_names = ["timfreq_bunch"],
        function = timfreqDev
    ),
    name = "bunch_log"
)
bunch_log.inputs.time_groups = c.groups
bunch_log.inputs.time_pool   = c.pooling
bunch_log.inputs.time_binary = c.binary

# Add regressors to Bunch
bunch_reg = pe.Node(
    Function(
        input_names = ["bunch", "confounds_path", "confounds_names"],
        output_names = ["design_bunch"],
        function = add_nuisance
    ),
    name = "bunch_reg"
)
bunch_reg.inputs.confounds_names = c.tapas_cols #TODO: read from physio.mat > hardcoding

# Unzip functional images (preprocessed BOLD)
unzip = pe.MapNode(
    Gunzip(),
    name = 'unzip',
    iterfield = ["in_file"] 
)

# --------- A. Generate design information - Specify the SPM model
# https://nipype.readthedocs.io/en/1.1.0/users/model_specification.html
modeler = pe.Node(
    model.SpecifySPMModel(
        concatenate_runs = c.concat,
        input_units  = 'secs',
        output_units = 'secs',
        high_pass_filter_cutoff = 128,
        parameter_source = "FSL"
    ),
    name = 'modeler'
)

# --------- B. Level1Design - Generate an SPM design matrix
designer = pe.Node(
    spm.Level1Design(
        bases = {'hrf': {'derivs': c.hrf_dervs}},
        timing_units = 'secs',
        volterra_expansion_order = (2 if c.volterra else 1)
    ),
    name = 'designer'
)

# --------- C. Estimate Model - Estimate the parameters of the model
estimator = pe.Node(
    spm.EstimateModel(estimation_method = {'Classical': 1}),
    name = 'estimator'
)

# --------- D. Contrastor -  Estimate contrasts
if c.contrast:
    contrastor = pe.Node(
        spm.EstimateContrast(contrasts = c.contrasts[c.jobName]),
        name = 'contrastor'
    )

# -------------------------------------------------------------------------------------------------
# 04. Connect the Nodes: Determine the Flow of Data
# -------------------------------------------------------------------------------------------------
timDev22 = Workflow(name = "l1_timDev")
timDev22.base_dir = str(c.workDir)

# Specify how the analysis iterates through the data
timDev22.connect([(infosource, infohandle, [
    ("subID", "subID"),
	("sesID", "sesID"),
	("acqID", "acqID")
    ])])

# Ñarse the logfiles into bunches
timDev22.connect([
    (infohandle, bunch_log, [("log_path", "time_log")]),
    (infohandle, bunch_reg, [("reg_path", "confounds_path")]) # reg_path: only exists when including BIOPAC
])

timDev22.connect([
    (bunch_log, bunch_reg, [("timfreq_bunch", "bunch")])
])

# Unzip bold files
timDev22.connect([(infohandle, unzip, [("bold_path", "in_file")])])

# Model specs
if c.smoothing is not None:
    timDev22.connect([
        (unzip, smoother, [("out_file", "in_files")]),
        (smoother, modeler, [("smoothed_files", "functional_runs")])
    ])
else:
    timDev22.connect([
        (unzip, modeler, [("out_file", "functional_runs")])
    ])

if c.artDetect:
    # Estimate motion outliers
    timDev22.connect([(infohandle, art_detect, [
            ("mask_path", "mask_file"),
            ("movpar_path", "realignment_parameters"),
            ("bold_path", "realigned_files")])])

    # Connect to modeler
    timDev22.connect([(art_detect, modeler, [("outlier_files", "outlier_files")])])
else:
    timDev22.connect([
        (infohandle, modeler, [
            ("out_path", "outlier_files"),
            ("conf_path", "realignment_parameters")])])

# Connect to modeler
timDev22.connect([
    (infohandle, modeler, [("TR", "time_repetition")]),
    (bunch_reg, modeler, [("design_bunch", "subject_info")])
])

# Design the matrix
timDev22.connect([
    (modeler, designer, [("session_info", "session_info")]),
    (infohandle, designer, [("TR", "interscan_interval")])])

# Estimate
timDev22.connect([(designer, estimator, [("spm_mat_file", "spm_mat_file")])])

# Contrast
timDev22.connect([
				(estimator, contrastor, [("spm_mat_file", "spm_mat_file")]),
				(estimator, contrastor, [("beta_images", "beta_images")]),
				(estimator, contrastor, [("residual_image", "residual_image")])])

# Save files
timDev22.connect([
    (estimator, datasink_T1w, [
        ('spm_mat_file', '1stLevel.@estimator_spm_mat'),
        ('residual_image', '1stLevel.@residuals'),
        ('beta_images', '1stLevel.@beta_images')
    ]),
    (contrastor, datasink_T1w, [
        ('spm_mat_file', '1stLevel.@contrastor_spm_mat'),
        ('spmT_images', '1stLevel.@T'),
        ('con_images', '1stLevel.@con')
    ])
])

# -------------------------------------------------------------------------------------------------
# 05. Visualize the Workflow
# -------------------------------------------------------------------------------------------------
timDev22.write_graph(graph2use='colored', format='png', simple_form=True)

# -------------------------------------------------------------------------------------------------
# 06. Run the Workflow
# -------------------------------------------------------------------------------------------------
res = timDev22.run()