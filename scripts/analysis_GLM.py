#! /usr/bin/env python
# Time-stamp: <2026-05-26 m.utrosa@bcbl.eu>
'''
Fixed effects fMRI model fitting
Concatenating time courses at corresponding locations across subjects

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
from nipype.interfaces.nipy.model import FitGLM
from nipype.interfaces.freesurfer import MRIConvert
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.io import DataSink
from nipype import Workflow, Function
import nipype.interfaces.spm as spm  # spm
import nipype.interfaces.matlab as mlab  # how to run matlab
import nipype.pipeline.engine as pe  # pypeline engine
# import nipype.algorithms.rapidart as ra  # artifact detection
import nipype.algorithms.modelgen as model  # model specification

# Import custom-made functions (scripts)
import grabber
from objects_timDev import grab_objects
from designs_timDev import timDev

# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")

# Set up project root, needed paths and folders
homePath = Path('/home/mutrosa/mutrosa/Documents/devLoc')
mriPath  = homePath / "data_MRI" / "derivatives" / "NORDIC-False" # path to preproc outputs
work_dir = homePath / "results" / "work" # for intermediate outputs
out_dir  = homePath / "results"
MNI      = homePath / "templates" / "tpl-MNI152NLin2009cAsym_res-01_T1w.nii.gz" # the same as in fMRIprep !

# Set up experimental procedure info
subjects = [5]
sessions = [2]
runs     = ['BLOCK1', 'BLOCK2', 'BLOCK3', 'BLOCK4']
anat_ses = 2

# Experimental design and MRI info
task = "timDev"

# Model parameters
smoothing  = None # Set the Gaussian filter width in mm, default is None

# T1w Datasink: create output folder for important outputs in T1w space
datasink_T1w = pe.Node(
    DataSink(
        base_directory = str(work_dir),
        container = str(out_dir)
    ),
    name = "datasink_T1w"
)

# MNI Datasink: create output folder for important outputs in MNI space
datasink_MNI = pe.Node(
    DataSink(base_directory = str(work_dir),
                             container = str(out_dir)),
    name = "datasink_MNI"
)

# Define a Node that extracts filepaths for all files required for the analysis
infohandle = pe.Node(
    Function(
        input_names  = ["subID", "sessions", "anatID", "runs", "homePath", "mriPath"],
        output_names = [
    "log_paths",
    "bold_paths", 
    "mask_paths", 
    "conf_paths",
    "out_paths", 
    "T1w_paths", 
    "T1w_toMNI_paths",
    "orig_to_boldref_paths",
    "boldref_to_T1w_paths", 
    "TRs"
    ],
        function = grab_objects),
name = "infohandle"
    )

infohandle.inputs.subID    = subjects[0]
infohandle.inputs.sessions = sessions
infohandle.inputs.anatID   = anat_ses
infohandle.inputs.runs     = runs
infohandle.inputs.homePath = str(homePath)
infohandle.inputs.mriPath  = str(mriPath)
# -------------------------------------------------------------------------------------------------
# 01. Specify 1st-level model parameters
# -------------------------------------------------------------------------------------------------
# Get the information about the experimental paradigm to create an SPM design matrix.
# Construct a list of objects (each object should contain data for all runs of that session)
# Create a Bunch object by parsing all event files of the
design_bunch = pe.Node(
    Function(
        input_names = ["logfilepaths"],
        output_names = ["design_info_list"],
        function = timDev
    ),
    name = "design_bunch"
)

# Unzip functional images (preprocessed BOLD)
unzip = pe.MapNode(
    Gunzip(),
    name = 'unzip',
    iterfield=['in_file']
)

# --------- A. Generate design information - specify the model
model_spec = pe.Node(
    interface = model.SpecifySPMModel(),
    name = "modelspec"
)
model_spec.inputs.concatenate_runs = True # treat runs as a single continuous series (fixed effects)!
model_spec.inputs.input_units = 'secs'
model_spec.inputs.output_units = 'secs'
model_spec.inputs.high_pass_filter_cutoff = 128 # High filter (default)

# --------- B. Fit the GLM model using nipy and ordinary least square method
model_estimate = pe.Node(
    interface = FitGLM(),
    name = "model_estimate"
    )
model_estimate.inputs.model = "spherical"
model_estimate.inputs.method = "ols"


# -------------------------------------------------------------------------------------------------
# 02. Connect the Nodes: Determine the Flow of Data
# -------------------------------------------------------------------------------------------------
timDev22 = Workflow(name = "level1")
timDev22.base_dir = str(work_dir)

timDev22.connect([

    # Generate lists for concatenation
    (infohandle, design_bunch, [("log_paths", "logfilepaths")]),

    # Generate lists of preprocesed data
    (infohandle, unzip, [("bold_paths", "in_file")]),

    # Model specs
    (design_bunch, model_spec, [("design_info_list", "subject_info")]),
    (infohandle, model_spec, [
            ("out_paths", "outlier_files"),
            ("conf_paths", "realignment_parameters"),
            ("TRs", "time_repetition") # CHECK: TR is a single float
            ])
    ])

if smoothing is not None:
    smooth = pe.Node(interface=spm.Smooth(), name="smooth")
    smooth.inputs.fwhm = smoothing
    timDev22.connect([
        (unzip, smooth, [("out_file", "in_files")]),
        (smooth, model_spec, [("smoothed_files", "functional_runs")])
    ])
else:
    timDev22.connect([
        (unzip, model_spec, [("out_file", "functional_runs")])
    ])

# Estimation
timDev22.connect([
    (infohandle, model_estimate, [("TRs", "TR")]),
    (model_spec, model_estimate, [('session_info', 'session_info')])
    ])
])
# -------------------------------------------------------------------------------------------------
# 03. Visualize the Workflow
# -------------------------------------------------------------------------------------------------
timDev22.write_graph(graph2use = 'colored', format = 'png', simple_form = True)

# -------------------------------------------------------------------------------------------------
# 04. Run the Workflow
# -------------------------------------------------------------------------------------------------
res = timDev22.run()
