#! /usr/bin/env python
# Time-stamp: <2026-06-04 m.utrosa@bcbl.eu>
'''
Fixed effects fMRI model fitting

Fixed-Effects Analysis
https://www.brainvoyager.com/bv/doc/UsersGuide/StatisticalAnalysis/FixedEffectsRandomEffectsMixedEffects.html
"The simple concatenation approach constitutes a fixed effects (FFX) analysis assessing observed activation
 effects with respect to the scan-to-scan measurement error, i.e. with respect to the precision with which 
 we can measure the fMRI signal. The source of variability used in a FFX analysis, thus, represents 
 within-subject variance."
- concatenation approach:  the number of analyzed data points is the sum of the data points from all sessions
- all N sessions have the same number of data points n, the total number of data points NT is NT = N x n.
- cannot be generalized to the population level (single-case study)

GLM in SPM
https://nipype.readthedocs.io/en/latest/users/examples/fmri_nipy_glm.html
'''

# Import python packages
from pathlib import Path
import pandas as pd
import numpy as np

# Import nipype stuff
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.io import DataSink
from nipype import Workflow, Function, IdentityInterface
import nipype.interfaces.spm as spm  # spm
import nipype.pipeline.engine as pe  # pypeline engine
import nipype.algorithms.modelgen as model  # model specification
# from nipype.algorithms.rapidart import ArtifactDetect # artifact detection

# Import custom-made functions (scripts)
import grabber
from objects import grab_objects
from designs import timfreqDev

# Set up project root, needed paths and folders
homePath = Path('/home/mutrosa/mutrosa/Documents/devLoc')
mriPath  = homePath / "data_MRI" / "derivatives" / "NORDIC-True" # path to preproc outputs
work_dir = homePath / "results" / "work" # for intermediate outputs
out_dir  = homePath / "results"

# Set up experimental procedure info
subjects =[5]
sessions = [2, 3, 4, 5, 6, 7]
acquisitions = ['BLOCK1', 'BLOCK2', 'BLOCK3', 'BLOCK4']
anat_ses = 2
space = "T1w"
pooling = True # If True absolute timing deviancy regressors, if False nominal.
hrf_dervs = [0, 0]
volterra = False

# Experimental design and MRI info
task = "timDev"

# Model parameters
smoothing  = 2 # Set the Gaussian filter width in mm. Default None (no smoothing).

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
    ('subID', subjects),
	('sesID', sessions),
	('acqID', acquisitions)
]

# T1w Datasink: create output folder for important outputs in T1w space
datasink_T1w = pe.Node(
    DataSink(
        base_directory = str(work_dir),
        container = str(out_dir)
    ),
    name = "datasink_T1w"
)

# Output substitutions: correct all Datasink output folder structures
substitutions = []
subjFolders = [('_acqID_%s_sesID_%s_subID_%s' % (acq, ses, sub),
				'sub-0%s/ses-0%s/acq-%s' % (sub, ses, acq))
               for acq in acquisitions
               for ses in sessions
               for sub in subjects]
substitutions.extend(subjFolders)
datasink_T1w.inputs.substitutions = substitutions
datasink_T1w.inputs.substitutions += [('beta_', f'beta_space-{space}_'),]

# Define a Node that extracts filepaths for all files required for the analysis
infohandle = pe.Node(
    Function(
        input_names  = [
            "subID",
            "sesID", 
            "anatID", 
            "homePath", 
            "mriPath",
            "space",
            "acqID", 
            "runs"
        ],
        output_names = [
            "log_path",
            "bold_path", 
            "mask_path", 
            "conf_path",
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
infohandle.inputs.anatID   = anat_ses
infohandle.inputs.homePath = str(homePath)
infohandle.inputs.mriPath  = str(mriPath)
infohandle.inputs.space    = space
infohandle.inputs.runs     = False

# -------------------------------------------------------------------------------------------------
# 02. Additional preprocessing node: smoothing and outlier detection
# -------------------------------------------------------------------------------------------------
if smoothing:
    smooth = pe.Node(interface=spm.Smooth(), name="smooth")
    smooth.inputs.fwhm = smoothing # TODO: Check whether this has to be a list?  

# Using intensity and motion parameters to infer parameters
# CHECK: which threshod to use? is this valid given that fMRIprep realigns data (does it?)?
# art_detect = pe.MapNode(
#     ArtifactDetect(),
#     name = "art_detect",
#     iterfield = ['in_file']
# )

# art_detect.inputs.realignment_parameters = 'functional.par'
# art_detect.inputs.parameter_source = 'FSL'
# art_detect.inputs.norm_threshold = 1
# art_detect.inputs.use_differences = [True, False]
# art_detect.inputs.zintensity_threshold = 3

# TODO: Detect outliers
# TODO: (art_detect, modeler, [("outlier_files", "outlier_files")]), # Optional: use fMRIprep/Nipype-computed artifacts]),
# (infohandle, art_detect, [
#     ("???", "realignment_parameters")
#     ("bold_paths", "realigned_files"),
#     ("mask_paths", "mask_type") # CHECK: are these correct masks?
#     ]),

# -------------------------------------------------------------------------------------------------
# 03. Specify 1st-level model parameters
# -------------------------------------------------------------------------------------------------
# Get the information about the experimental paradigm to create an SPM design matrix.
# Construct a list of objects (each object should contain data for all runs of that session)
# Create a Bunch object by parsing all event files: timDev & freqDev are separate Bunch objects.
design_bunch = pe.Node(
    Function(
        input_names = ["time_log", "time_pool"],
        output_names = ["timfreq_bunch"],
        function = timfreqDev
    ),
    name = "design_bunch"
)
design_bunch.inputs.time_pool = pooling

# Unzip functional images (preprocessed BOLD)
unzip = pe.MapNode(
    Gunzip(),
    name = 'unzip',
    iterfield = ["in_file"] 
)

# --------- A. Generate design information - Specify the SPM model
modeler = pe.Node(
    model.SpecifySPMModel(
        concatenate_runs = False, # treat runs as a single continuous series (fixed effects)!
        input_units  = 'secs',
        output_units = 'secs',
        high_pass_filter_cutoff = 128),
    name = 'modeler'
)

# --------- B. Level1Design - Generate an SPM design matrix
designer = pe.Node(
    spm.Level1Design(
        bases = {'hrf': {'derivs': hrf_dervs}},
        timing_units = 'secs',
        volterra_expansion_order = (2 if volterra else 1)
    ),
    name = 'designer'
)

# --------- C. Estimate Model - Estimate the parameters of the model.
estimator = pe.Node(
    spm.EstimateModel(estimation_method = {'Classical': 1}),
    name = 'estimator'
)

# -------------------------------------------------------------------------------------------------
# 02. Connect the Nodes: Determine the Flow of Data
# -------------------------------------------------------------------------------------------------
timDev22 = Workflow(name = "level1")
timDev22.base_dir = str(work_dir)
timDev22.connect([(infosource, infohandle, [
    ("subID", "subID"),
	("sesID", "sesID"),
	("acqID", "acqID")
    ])
])
timDev22.connect([

    # Generate lists for concatenation
    (infohandle, design_bunch, [("log_path", "time_log")]),

    # Generate lists of preprocesed data
    (infohandle, unzip, [("bold_path", "in_file")])
])

if smoothing is not None:
    timDev22.connect([
        (unzip, smooth, [("out_file", "in_files")]),
        (smooth, modeler, [("smoothed_files", "functional_runs")])
    ])
else:
    timDev22.connect([
        (unzip, modeler, [("out_file", "functional_runs")])
    ])

timDev22.connect([

    # Model specs
    (infohandle, modeler, [
            ("out_path", "outlier_files"),
            ("conf_path", "realignment_parameters"),
            ("TR", "time_repetition")
    ]),
    (design_bunch, modeler, [("timfreq_bunch", "subject_info")])
])

timDev22.connect([
    (modeler, designer, [("session_info", "session_info")]),
    (infohandle, designer, [("TR", "interscan_interval")])
])

timDev22.connect([
    (designer, estimator, [("spm_mat_file", "spm_mat_file")]),
    (estimator, datasink_T1w, [
        ('spm_mat_file', '1stLevel.@spm_mat'),
        ('residual_image', '1stLevel.@residuals'),
        ('beta_images', '1stLevel.@beta_images')
    ]),
])

# -------------------------------------------------------------------------------------------------
# 03. Visualize the Workflow
# -------------------------------------------------------------------------------------------------
timDev22.write_graph(graph2use = 'colored', format = 'png', simple_form = True)

# -------------------------------------------------------------------------------------------------
# 04. Run the Workflow
# -------------------------------------------------------------------------------------------------
res = timDev22.run()