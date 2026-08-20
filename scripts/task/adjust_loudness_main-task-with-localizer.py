#! /usr/bin/env python
# Time-stamp: <2026-04-23 m.utrosa@bcbl.eu>
'''
Script to test loudness of sounds:
- localizer: naturalistic sounds
- main task: complex harmonic sounds
'''
# 00. PREPARATION ----------------------------------------------------------------------------------
import ast
import random
import pandas as pd
import sounddevice as sd
from pathlib import Path
from expyriment import design, control, stimuli, misc, io

# Import the external sequence generation file for the main task
import create_soundtrack_soundgen as sg

# 01. PARAMETERS -------------------------------------------------------------------------
sesID = 27
control.set_develop_mode(on=True)
main_task_on = True
params = {

	# Local setup
    "PROJECT_ROOT" : "C:/Users/Experimental User/Desktop/SUBCORT_HIGHRES/",
    "AUDIO_ROOT" : "C:/Users/Experimental User/Desktop/SUBCORT_HIGHRES/stimuli",
    "AUDIOFILE_REGEX" : "**/*.wav",
    
    # Visual
    "CANVAS_SIZE"             : (1024, 768), # MRI monitor resolution.
    "FIXATION_CROSS_SIZE"     : (40, 40),
    "FIXATION_CROSS_POSITION" : (0, 0),
    "FIXATION_CROSS_WIDTH"    : 4,

    # Colors in RGB
    "BLACK"   : (0, 0, 0),       # screen background
    "WHITE"   : (255, 255, 255), # fixation cross

    # Audio for main task
    "TONE_LOUDNESS"   : 85,     # dB SPL
    "TONE_DURATION"   : 50,     # msec
    "NUM_HARMONICS"   : 10,     # Number of harmonics
    "HARMONIC_FACTOR" : 0.8,    # Harmonic amplitude decay factor
    "MAX_AMPLITUDE"   : 1.14,   # Defined through a simulation
    "SAMPLE_RATE"     : 48000,  # Hz
    "TAU"             : 5,      # Ramping window in msec

	# Audio for localizer
	"SOUND_DURATION"   : 1000, # msec
	"SOUNDS_PER_TRIAL" : 2,
}

# 02. PREPARATION ----------------------------------------------------------------------------------
# Load audio stimuli for localizer
audio_root    = Path(params["AUDIO_ROOT"])
wav_filepaths = list(Path(params["AUDIO_ROOT"]).glob(params["AUDIOFILE_REGEX"]))

# Rule: sounds include "s3" in filename & silences "null".
filenames_sounds = [file_1 for file_1 in wav_filepaths if "s3" in str(file_1)]
filename_null    = [file_2 for file_2 in wav_filepaths if "null" in str(file_2)][0]

# Shuffle the sounds
random.shuffle(filenames_sounds)

# Load the trial parameters from csv for main task
homePath  = Path(params["PROJECT_ROOT"])
paramPath = homePath / f"ses-{sesID:003d}_exp_parameter_combo.csv"
df        = pd.read_csv(paramPath)
no_blocks = len(df["block_no"].unique())

# Ensure that the trials are ordered by block & trial IDs
df.sort_values(by=["block_no", "trial_no"], inplace=True)

# Ensure correct data types in columns with lists as row values.
list_cols = ["freq_dev", "freq_dev_type", "freq_loc", "freq_diff", "freq_diff_abs"]
for col in list_cols:
    df[col] = df[col].apply(
    lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )

# 03. INITIALIZE THE EXPERIMENT ----------------------------------------------------------
control.defaults.window_size=params["CANVAS_SIZE"]
control.defaults.window_no_frame=True
control.defaults.opengl=2
control.defaults.display_resolution=params["CANVAS_SIZE"]
exp  = design.Experiment(name = "loudness")
control.initialize(exp)

# 04. CREATE & PRELOAD THE STIMULI -----------------------------------------------------------------
# Creating
keyboard = io.Keyboard()
canvas = stimuli.Canvas(size=params["CANVAS_SIZE"], colour=params["BLACK"])
fix_cross = stimuli.FixCross(size=params["FIXATION_CROSS_SIZE"], position=params["FIXATION_CROSS_POSITION"], 
								  line_width=params["FIXATION_CROSS_WIDTH"], colour=params["WHITE"])
# Preloading
fix_cross.preload(); fix_cross.plot(canvas); canvas.preload();

# Get sounds for localizer for directory
sounds = {filename: stimuli.Audio(str(filename)) for filename in filenames_sounds}
for s in sounds.values():
    s.preload()

# Get sounds for main task: initialize soung generation (SoundGen) class
sound_gen = sg.SoundGen(params["SAMPLE_RATE"], params["TAU"])

# 05. RUN THE EXPERIMENT ---------------------------------------------------------------------------
# Start the loudness adjustment.
control.start(skip_ready_screen=True) # Start the experiment without the ready screen and wait for trigger from the MRI.
canvas.present()

### LOCALIZER SOUNDS
while True:
	key, rt = keyboard.wait(keys = [misc.constants.K_g, misc.constants.K_e])

	# Play sounds when 'g' is pressed (GO).
	if key == 103: # ASCII code
		for i in range(params['SOUNDS_PER_TRIAL']):
			i_random = random.choice(list(sounds.items()))
			i_sound = i_random[1]
			i_sound.play()
			exp.clock.wait(params["SOUND_DURATION"])

	# End loudness check when 'e' is pressed (END).
	if key == 101: # ASCII code
		break

### MAIN TASK SOUNDS
for block in range(no_blocks):

    # Correct for zero-indexing
    block_idx = block + 1

    # Select the part of the dataframe relevant for the current block
    df_block = df[df["block_no"] == block_idx]

    # Mark the start of the functional sequence for the main task
    task_start_time = exp.clock.time

    # Play all tone sequences: trial by trial
    block_start_time = exp.clock.time - task_start_time
    for soundtrack, ITI, freq_dev_no, trial_log in sound_gen.generate_soundtrack(df_block, block_start_time, params["MAX_AMPLITUDE"], params["NUM_HARMONICS"],  params["TONE_DURATION"],  params["HARMONIC_FACTOR"], params["TONE_LOUDNESS"]):

        key, rt = keyboard.wait(keys = [misc.constants.K_g, misc.constants.K_e])

        # Play sounds when 'g' is pressed (GO).
        if key == 103: # ASCII code
            sd.play(soundtrack, samplerate = params["SAMPLE_RATE"])
            sd.wait()
        
        # End loudness check when 'e' is pressed (END).
        if key == 101: # ASCII code
            # Finishing
            control.end()
