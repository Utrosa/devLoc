#! /usr/bin/env python
# Time-stamp: <27-04-2026, m.utrosa@bcbl.eu>
'''
Main experimental script
- runs localizer
- runs frequency counting task with timing deviant tone sequences
- logs events in BIDS format

Requires Expyriment version 1.0.0
Local setup: All localizer sounds must be in "Stimuli" folder. Sounds need "s3" prefix.
'''

# 00. PREPARATION ------------------------------------------------------------------------
import ast
import random
import numpy as np
import pandas as pd
import sounddevice as sd
from pathlib import Path
from datetime import datetime 
from expyriment import design, control, stimuli, misc, io

# Import the external sequence generation file for the main task
import create_soundtrack_soundgen as sg

# Specify BIDS-formatted EventFiles for localizer
## onset [sec], duration [sec], stim_file [wav],
## key [chr(ASCII)], RT [sec]
log_loc_format_fStr  = "{0:.3f}\t{1:.3f}\t{2}\t{3}\t{4:.3f}\n"
log_loc_format_NaNs  = "{0:.3f}\t{1:.3f}\t{2}\t{3}\t{4}\n"
log_task_format_fStr = "{1}\t{2}\t{3}\t{4}\n"

# 01. PARAMETERS -------------------------------------------------------------------------
idx  = input("Enter the index of exp_parameter_combo csv:")
sesh = input("Enter the session number:")
comboID = int(idx); sesID = int(sesh)
control.set_develop_mode(on=False)
localizer_on = True
main_task_on = True
params = {

    # Directories
    "PROJECT_ROOT" : "C:/Users/Experimental User/Desktop/SUBCORT_HIGHRES",
    "AUDIO_ROOT"   : "C:/Users/Experimental User/Desktop/SUBCORT_HIGHRES/stimuli",
    "AUDIOFILE_REGEX" : "**/*.wav",
    
    # Experiment structure
    "WAIT_TIME" : 3000, # msec (ensuring reading at the start & end of the experiment)
    
    # Localizer structure
    "LOC_REP"             : 1,  # number of repetitions
    "LOC_TRIALS"          : 10, # the number of equally long sound and silence pairs
    "SOUNDS_PER_SEQUENCE" : 30, # determines the length of trials; each sound is 1 sec

    # Visual
    "CANVAS_SIZE"             : (1280, 800), # MRI monitor resolution.
    "FIXATION_CROSS_SIZE"     : (40, 40),
    "FIXATION_CROSS_POSITION" : (0, 0),
    "FIXATION_CROSS_WIDTH"    : 5,
    "HEADING_SIZE"            : 30, 
    "TEXT_SIZE"               : 25, 
    
    # Colors in RGB
    "BLACK"   : (0, 0, 0),       # screen background
    "WHITE"   : (255, 255, 255), # fixation cross
    "CORRECT" : (70, 255, 255),  # cyan; color contrast on black is super (17.01)
    "WRONG"   : (255, 234, 0),   # yellow; color contrast on black is super (17.02)
    
    # Audio for task
    "TONE_LOUDNESS"   : 85,     # dB SPL
    "TONE_DURATION"   : 50,     # msec
    "NUM_HARMONICS"   : 10,     # Number of harmonics
    "HARMONIC_FACTOR" : 0.8,    # Harmonic amplitude decay factor
    "MAX_AMPLITUDE"   : 1.14,   # Defined through a simulation
    "SAMPLE_RATE"     : 48000,  # Hz
    "TAU"             : 5,      # Ramping window in msec
    
    # Sound stimuli in localizer
    "SOUND_STRATA"     : 84,   # the total amount of available sounds
    "SOUND_DURATION"   : 1000, # msec
    "SOUND_REP_PROB"   : .05,  # low enough to be stimulating/challenging

    # In the MRI, the rainbow response pad is 1234, while the gun one abcd.
    "DETECTION_SYMBOL"   : [misc.constants.K_1, misc.constants.K_2, misc.constants.K_3, misc.constants.K_4],
    "DETECTION_ASCII"    : ['49', '50', '51', '52'],

    # Text (for instructions and goodbye message)
    "INTRO_HEADING"   : f"Welcome to the ",
    "INTRO_TEXT_TASK" : f"When the experiment starts, you will see a white cross in the center of the screen.\n"
                        "Please look at the cross for the entire duration of the task.\n\n"
                        "You will hear a series of sounds. Each series has 7 sounds.\n"
                        "Your task is to count how many sounds have a different pitch compared to the rest.\n\n "
                        "There may be 0, 1, 2 or 3 sounds with a different pitch.\n"
                        "After you hear the last sound of the series, press 0, 1, 2 or 3.\n"
                        "Try to be as fast and as accurate as possible.\n\n"
                        "If you counted well and were fast enough, the cross will change color to cyan.\n"
                        "If you made a mistake or were too slow, the cross will change color to yellow.\n\n"
                        "Sometimes, there will be no sound for a while.\n"
                        "During these periods, please just relax and stay still.\n\n"
                        "We will repeat this task four times.\n"
                        "You will have short breaks in between.\n\n\n"
                        "It is very important that you do NOT move your arms, legs, or head until the break.\n\n"
                        "Thank you and good luck!\n\n\n"
                        "Whenever you are ready, press any button.",
    "INTRO_TEXT_LOC"  :  f"When the experiment starts, you will see a white cross in the center of the screen.\n"
                        "Please look at the cross for the entire duration of the task.\n\n"
                        "You will hear a series of short sounds (animals, speech, or tools).\n"
                        "Sometimes, the same sound will play twice in a row.\n"
                        "When you hear a sound repeat, press any button as fast as you can.\n\n"
                        "If you you were correct and fast enough, the cross will change color to cyan.\n"
                        "If you made a mistake or were too slow, the cross will change color to yellow.\n\n"
                        "Each sound series will be separated by a longer silent period.\n"
                        "During these periods, please just relax and stay still.\n\n"
                        "It is very important that you do NOT move your arms, legs, or head until the break.\n\n"
                        "Thank you so much and good luck!\n\n\n"
                        "Whenever you are ready to start, press any button.",

    "REST_HEADING" : "BREAK TIME",
    "REST_TEXT"    : "\n\nPlease take a moment to rest and move your body as needed.\n\n"
                     "Whenever you are ready to continue, press any button.",
    
    "MRI_HEADING"  : "SCANNER CALIBRATION",
    "MRI_TEXT"     : "Please remain still for a few moments.\n\n"
                     "Remember NOT to move your body or head until the break.\n\n"
                     "Thank you!",

    "TASK_END_HEADING" : "You have reached the end of the ", 
    "TASK_END_TEXT"    : "Thank you so much!\n\n",

    "END_HEADING"      : "The End of the Experiment", 
    "END_TEXT"         : "Thank you so much for your participation!\n\n",
}

# 2. FUNCTIONS --------------------------------------------------------------------------
# Functions needed for the localizers.
def create_soundtrack(sound_strata, sequence_len, rep_prob, sequence_no):
    '''
    Generates sounds sequences with the following constraints:
    (a) Each sound has a probability of repetition defined by rep_prob.
    (b) The last two and first three sounds in the experiment are always unique.
    (c) Sounds are distributed randomly across the experiment.
    (d) All sounds are used fairly, avoiding selection bias.

    Parameters:
    - sound_strata: A list of all loaded sounds.
    - sequence_len: The number of sounds per sequence.
    - rep_prob: The probability of a sound repeating within the experiment.
    - sequence_no: The total number of sequences.

    Returns:
    - sequences: A list of lists (sound sequences) with the above constraints.
    '''

    # Calculate the number of unique and repeated sounds across the sequences.
    total_sounds = sequence_len * sequence_no
    n_reps_float = rep_prob * total_sounds
    n_reps = int(np.floor(n_reps_float) + (np.random.rand() < (n_reps_float - np.floor(n_reps_float))))
    n_norep_max = total_sounds - n_reps
    if n_reps > total_sounds - 5:
        raise ValueError("Probability of repetitions is too high.")

    # If the number of sounds needed in total is larger then the number of available sounds,
    # prolong the sound strata by looping through it multiple times.
    n_loops = int(np.ceil(total_sounds / len(sound_strata)))

    # Create a list of non-repeated sounds. A repetition is when two sounds repeat consecutively.
    resample = True
    while resample:
        all_strata_sounds = [sound_strata[u] for _ in range(n_loops) for u in np.random.permutation(len(sound_strata))]
        
        # Check for sequential repetitions in the generated sequence: no previous and current sound are the same.
        resample = any([all_strata_sounds[u] == all_strata_sounds[u - 1] for u in range(1, len(all_strata_sounds))])
    
    all_unique_sounds = all_strata_sounds[:n_norep_max] # Shorten to the desired number of unique sounds.
    all_sounds = all_strata_sounds[:total_sounds] # Shorten to the desired length of the experiment.

    # Determine the indices of both sequence clashes (boundaries) and necessarily unique sounds.
    all_idx = list(range(0,len(all_sounds)))
    unique_idx = [0, 1, 2] + all_idx[-2:] # No repetitions for the last 2 and first 3 sounds in the experiment.
    boundary_idx = [[all_idx[sequence_len * n - 1]] + [all_idx[sequence_len * n]] for n in range(sequence_no)]
    boundary_idx = boundary_idx[1:] # The first pair is irrelevant (the indices of the first and last sound in the exp.)
    
    # Flatten the boundary indices list and add it to the unique_idx list.
    unique_idx = np.sort(unique_idx + [item for sublist in boundary_idx for item in sublist])

    # Determine the sound repetition indices (two values per repetition).
    resample = True
    while resample:

        # Randomly sample repetition indices, ensuring no overlap with unique_idx.
        repeat_idx = [(idx, idx + 1) for idx in np.sort(np.random.permutation(range(3, total_sounds - 1))[:n_reps])]
        flattened_repeat_idx = [x for idx_tuple in repeat_idx for x in idx_tuple]

        # Resample if idx duplicates or if any idx is in both repeat_idx and unique_idx.
        resample = any([ix0 == ix1 for ix0, ix1 in zip(flattened_repeat_idx[1:], flattened_repeat_idx[:-1])])
        resample = resample or bool(set(flattened_repeat_idx) & set(unique_idx))

        assert len(repeat_idx) == n_reps # This should be unnecessary if the resampling logic is correct.

    # Create a sequence of all sounds presented in the experiment.
    for idx1, idx2 in repeat_idx:
        repeated_sound   = all_sounds[idx1]
        all_sounds[idx2] = repeated_sound

    # Separate the all_sounds sequence into separate sequences.
    sequences = [all_sounds[sequence_len * n : sequence_len * (n + 1)] for n in range(sequence_no)]

    return sequences

def compute_durations(pars, clock_start, clock_end, verbose=False):
    '''
    Computes the predicted and actual durations of sound or silence parts in trials.

    Parameters:
    - pars: A dictionary, containing all user-defined parameters.
    - clock_start: Start time given by the Expyriment's clock.
    - clock_end: Start time given by the Expyriment's clock.
    - verbose: Prints the computed durations in the terminal when True.

    Returns:
    - dur_predicted: The predicted trial duration, which is calculated from the user-defined audio parameters.
    - dur_actual: The actual trial duration, which is calculated from times given by Expyriment's clock.
    '''
    dur_predicted = pars['SOUNDS_PER_SEQUENCE'] * pars['SOUND_DURATION']
    dur_actual    = clock_end - clock_start

    if verbose:
        print(f'Predicted duration: {dur_predicted} msec. Actual duration: {dur_actual} msec.')

    return (dur_predicted, dur_actual)

def give_feedback(current_sound, position_in_sequence, response, good, bad):
    '''
    Determines the response category for a given response to a sound display.
    gives feedback by changing the color of the fixation cross.
    
    Parameters:
    - current_sound: The name of the sound that is currently being played (str).
    - position_in_sequence: The index of the sound in the current sound sequence (int).
    - response: The identity of the key pressed in ASCII code.
    - good: A class implementing fixation cross for correct responses (cyan).
    - bad: A class implementing fixation cross for wrong responses (yellow).

    Returns:
    - perfo_code: The response category (str).
    - run_performance: Dictionary tracking performance across runs (dict).
    - current_sound_key: The current stimuli name as obtained outside the funtion (str).
    - feeedback_status: Indicates whether feedback is given or not (boolean).
    '''
    # sounds_in_sequence is an empty dictionary. 
    # reversed_strata represents the reversed sound_strata dictionary. 
    # Dict keys are sounds names (e.g.: s3_animal_1_ramp10.wav) 
    # Dict values are expyriment audio object IDs, (e.g.: <expyriment.stimuli._audio.Audio object at 0x0000018DED268E08>).
    global sounds_in_sequence, reversed_strata, run_performance

    # Define current and previous sounds with names. 
    # Each sound presented in the experiment has a unique ID.
    # A sound repetition will have the same name but different ID.
    prev_sound_key    = sounds_in_sequence[position_in_sequence - 1] if position_in_sequence > 0 else None
    current_sound_key = reversed_strata.get(current_sound)

    # Default values
    perfo_code, feedback_status = None, None

    # Feedback structure
    if current_sound_key == prev_sound_key:
        if response is not None:
            good.present(clear = False, log_event_tag = True); run_performance["H"] += 1; perfo_code = "HIT"; feedback_status = True
        elif response is None:
            bad.present(clear = False, log_event_tag = True); run_performance["M"] += 1; perfo_code = "MISS"; feedback_status = True
    
    else:
        if response is None:
            run_performance["CR"] += 1; perfo_code = "CORR_REJECTION"; feedback_status = False
        elif response is not None:
            bad.present(clear = False, log_event_tag = True); run_performance["FA"] += 1; perfo_code = "FALSE_ALARM"; feedback_status = True

    return perfo_code, run_performance, current_sound_key, feedback_status

def play_sounds(sequence, sound_duration, exp, canvas, fixation, correct, wrong, keyboard, response_keys, log_events_sound):
    '''
    Plays one sound sequence, where sounds are presented one after the other.

    Parameters:
    - sequence: A single sequence of sounds (list) from the create_soundtrack() output.
    - sound_duration: The duration of individual sounds. All sounds need to have the same duration.
    - exp: A class implementing a basic experiment in expyriment.
    - canvas: A class implementing a canvas stimulus in expyriment.
    - fixation: A class implementing a general fixation cross in expyriment (white).
    - correct: A class implementing fixation cross for correct responses (cyan).
    - wrong: A class implementing fixation cross for incorrect responses (yellow).  
    - keyboard: A class implementing a keyboard input in expyriment.
    - response_keys: A list of accepted keys for responding (depends on the chosen response pad).
    - log_events_sound: A file containing information about events, following BIDS specification.

    Returns:
    - sounds_start: Time (float), indicating the start of the sound sequence.
    - sounds_end: Time (float), indicting the end of the sound sequence.
    '''
    global sounds_in_sequence, run_start_time, run_performance
    feedback_shown, fs = None, None
    sounds_start = exp.clock.time

    # Sound ID refers to the specific sound in the experiment
    # Sound repetitions have different sound IDs
    for count, sound_ID in enumerate(sequence):

        # Check if quit key is pressed
        keyboard.check(keys=[misc.constants.K_y])

        # Refresh the screen when feedback was given
        if feedback_shown is not None and count - feedback_shown == 2:
            canvas.present()

        # Stamp audio start time and play
        audio_start = exp.clock.time
        sound_ID.play(maxtime=sound_duration, log_event_tag=True)

        key_ASCII_audio, key_log_entry = None, None

        # Check for key presses while the sound is playing
        while sound_ID.is_playing and sound_ID.time < audio_start + sound_duration:

            # Non-blocking function to check for pressed keys
            keys = keyboard.read_out_buffered_keys()

            # ------ Key-press trials ------ 
            if keys and keys[0] != 115: # 115 is "s" from scanner sync box
                press_time = exp.clock.time
                key_ASCII_audio = keys[-1] # If multiple, we take the last key

                # Show feedback
                perf_code, run_performance, cs, fs = give_feedback(sound_ID,
                                                                   count,
                                                                   key_ASCII_audio,
                                                                   correct,
                                                                   wrong)
                
                # Logging
                if key_ASCII_audio is not None:
                    key_log_entry = {
                        "onset":      np.abs(audio_start - run_start_time) / 1000,
                        "duration":   sound_duration / 1000, # dummy duration
                        "stim_file":  cs.name,
                        "response":   chr(key_ASCII_audio),
                        "RT":         np.abs(press_time - audio_start) / 1000
                    }

        # After the sound has stopped playing
        audio_end = exp.clock.time

        # Correct sound duration for key trials
        if key_log_entry is not None:
            key_log_entry["duration"] = np.abs(audio_end - audio_start) / 1000
            log_events_sound.write(log_loc_format_fStr.format(
                                                key_log_entry["onset"],
                                                key_log_entry["duration"],
                                                key_log_entry["stim_file"],
                                                key_log_entry["response"],
                                                key_log_entry["RT"]
                                    ))

        # ------ No-key trials ------   
        if not key_ASCII_audio:
            
            # Show feedback
            perf_code, run_performance, cs, fs = give_feedback(sound_ID,
                                                               count,
                                                               key_ASCII_audio,
                                                               correct,
                                                               wrong)
            
            # Logging
            log_events_sound.write(log_loc_format_NaNs.format(
                np.abs(run_start_time - audio_start) / 1000, # onset
                np.abs(audio_end - audio_start) / 1000,      # duration
                cs.name,                                     # stim_file
                np.nan,                                      # response
                np.nan                                       # RT
            ))
        
        # Update feedback tracking and sequence record
        if fs:
            feedback_shown = count

        # Update sounds in sequence for correct feedback
        sounds_in_sequence.append(cs)

    sounds_end = exp.clock.time
    return sounds_start, sounds_end

def play_silence(null_sound, sound_duration, exp, null_number, keyboard, response_keys, log_events_null):
    '''
    Plays one silent sequence, where a single null sound is presented repetitively.

    Parameters:
    - null_sound: A single preloaded null sound.
    - sound_duration: The duration of the null sound.
    - exp: A class implementing a basic experiment in expyriment.
    - null_number: The number of null sound repetitions.
    - keyboard: A class implementing a keyboard input in expyriment.
    - response_keys: A list of accepted keys for responding (depends on the chosen response pad).
    - log_events_null: A file containing information about events, following BIDS specification.

    Returns:
    - silence_start: Time (float), indicating the start of the silence.
    - silence_end: Time (float), indicting the end of the silence.
    '''
    global run_start_time
    silence_start = exp.clock.time

    for null_event in range(null_number):

        # Check if quit key is pressed
        keyboard.check(keys=[misc.constants.K_y])

        # Refresh screen (needed if quit key is pressed)
        canvas.present()

        # Stamp audio start time and play
        null_start = exp.clock.time
        null_sound.play(maxtime=sound_duration, log_event_tag=True)

        key_ASCII_silence, key_log_entry = None, None

        # Check for key presses while the sound is playing
        while null_sound.is_playing and null_sound.time < null_start + sound_duration:

            # Non-blocking function to check for pressed keys
            keys = keyboard.read_out_buffered_keys()

            # ------ Key-press trials ------ 
            if keys and keys[0] != 115: # 115 is "s" from scanner sync box
                press_time = exp.clock.time
                key_ASCII_silence = keys[-1] # If multiple, we take the last key

                # Logging
                if key_ASCII_silence is not None:
                    key_log_entry = {
                        "onset": np.abs(null_start - run_start_time) / 1000,
                        "duration": sound_duration / 1000, # dummy duration
                        "stim_file": "null_event.wav",
                        "response": chr(key_ASCII_silence),  
                        "RT": np.abs(press_time - null_start) / 1000
                    }

        # After the sound has stopped playing
        null_end = exp.clock.time

        # Correct the sound duration for key trials
        if key_log_entry is not None:
            key_log_entry["duration"] = np.abs(null_end - null_start) / 1000
            log_events_null.write(log_loc_format_fStr.format(
                                                key_log_entry["onset"],
                                                key_log_entry["duration"],
                                                key_log_entry["stim_file"],
                                                key_log_entry["response"],
                                                key_log_entry["RT"]
                                    ))

        # ------ No-key trials ------ 
        if not key_ASCII_silence:
            log_events_null.write(log_loc_format_NaNs.format(
                np.abs(run_start_time - null_start) / 1000, # onset
                np.abs(null_end - null_start) / 1000,       # duration
                "null_event.wav",                           # stim_file
                np.nan,                                     # response
                np.nan                                      # RT
            ))
    
    silence_end = exp.clock.time
    return silence_start, silence_end

# Functions needed for the main task.
def give_feedback_freqCount(current_event, response, good, bad):
    '''
    Determines the response category for a given response in the frequency counting task.
    Gives feedback by changing the color of the fixation cross.
    
    Parameters:
    - current_event: The number of frequency deviants in the current trial.
    - response: The identity of the key pressed in ASCII code.
    - good: A class implementing fixation cross for correct responses.
    - bad: A class implementing fixation cross for wrong responses.

    Returns:
    - perfo_code: The response category (str).
    - trial_performance: Dictionary tracking performance across trials (dict).
    '''
    perfo_code = None

    # ----- Silent trials
    if np.isnan(current_event):
        current_event = np.nan_to_num(current_event)
    
    # Correct the response for MRI buttons: 1 is 0 devs, ... etc.
    response = response - 1

    # Ensure that you are comparing variables of the same type (str)
    current_event_str = str((int(current_event)))
    response_str = chr(response)

    # ----- Correct count -----
    if current_event_str == response_str:
        good.present(clear = False, log_event_tag = True)
        trial_performance["CORRECT"] += 1
        perfo_code = "CORRECT"
    
    # ----- Incorrect count ----- 
    else:
        bad.present(clear = False, log_event_tag = True)
        trial_performance["INCORRECT"] += 1
        perfo_code = "INCORRECT"

    return perfo_code, trial_performance

# 03. LOAD STIMULI -----------------------------------------------------------------------
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
paramPath = homePath / f"ses-{comboID:03d}_exp_parameter_combo.csv"
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

# 04. INITIALIZE THE EXPERIMENT ----------------------------------------------------------
# Settings to ensure that the window size fits the stimuli computer & projector in the MRI
control.defaults.opengl = 2 # OpenGL (vsync / blocking); default
control.defaults.display_resolution = params["CANVAS_SIZE"] # default is None
control.defaults.window_size = params["CANVAS_SIZE"] # relevant when running in window mode
control.defaults.window_no_frame = True              # relevant when running in window mode

# Initialize the experiment with a name following BIDS specification
exp  = design.Experiment(name = "devLoc")
control.initialize(exp)

# 05. CREATE & PRELOAD THE STIMULI -------------------------------------------------------
# Creating stimuli.
keyboard = io.Keyboard()
clock    = misc.Clock()

canvas        = stimuli.Canvas(size = params["CANVAS_SIZE"], colour = params["BLACK"])
blank_canvas  = stimuli.Canvas(size = params["CANVAS_SIZE"], colour = params["BLACK"])

scanner_text  = stimuli.TextScreen(
    params["MRI_HEADING"],
    params["MRI_TEXT"],
    heading_size=params["HEADING_SIZE"],
    heading_colour=params["WHITE"],
    text_size=params["TEXT_SIZE"],
    text_colour=params["WHITE"]
    )
instructions_loc  = stimuli.TextScreen(
    params["INTRO_HEADING"] + "repetition detection experiment!",
    params["INTRO_TEXT_LOC"],
    heading_size=params["HEADING_SIZE"],
    heading_colour=params["WHITE"],
    text_size=params["TEXT_SIZE"],
    text_colour=params["WHITE"]
    )
instructions_task = stimuli.TextScreen(
    params["INTRO_HEADING"] + "counting experiment!",
    params["INTRO_TEXT_TASK"],
    heading_size=params["HEADING_SIZE"],
    heading_colour=params["WHITE"],
    text_size=params["TEXT_SIZE"],
    text_colour=params["WHITE"]
    )
goodbye_exp_message = stimuli.TextScreen(
    params["END_HEADING"],
    params["END_TEXT"],
    heading_size=params["HEADING_SIZE"],
    heading_colour=params["WHITE"],
    text_size=params["TEXT_SIZE"],
    text_colour=params["WHITE"]
    )

fix_cross = stimuli.FixCross(size = params["FIXATION_CROSS_SIZE"], position = params["FIXATION_CROSS_POSITION"], line_width = params["FIXATION_CROSS_WIDTH"], colour = params["WHITE"])
wrong     = stimuli.FixCross(size = params["FIXATION_CROSS_SIZE"], position = params["FIXATION_CROSS_POSITION"], line_width = params["FIXATION_CROSS_WIDTH"], colour = params["WRONG"])
correct   = stimuli.FixCross(size = params["FIXATION_CROSS_SIZE"], position = params["FIXATION_CROSS_POSITION"], line_width = params["FIXATION_CROSS_WIDTH"], colour = params["CORRECT"])

# Get sounds for localizer for directory
silence = stimuli.Audio(str(filename_null))
sounds  = {filename: stimuli.Audio(str(filename)) for filename in filenames_sounds}

# Get sounds for main task: initialize soung generation (SoundGen) class
sound_gen = sg.SoundGen(params["SAMPLE_RATE"], params["TAU"])

# Preload to ensure fast stimuli presentation.
blank_canvas.preload(); scanner_text.preload()
instructions_loc.preload(); instructions_task.preload()
goodbye_exp_message.preload(); wrong.preload(); correct.preload()
silence.preload()
for s in sounds.values():
    s.preload()

# Plot fixation cross on canvas
fix_cross.preload(); fix_cross.plot(canvas)
canvas.preload()

# Create preloaded soundtrack for the localizer
## Pair sound ID names (for expyriment presentation) and filenames (for data storage).
sound_strata = random.sample(list(sounds.items()), params["SOUND_STRATA"])
sound_strata = dict(sound_strata)
reversed_strata = {value: key for key, value in sound_strata.items()} # Works because values are unique!
soundtrack = create_soundtrack(
    sound_strata = list(sound_strata.values()),
    sequence_len = params["SOUNDS_PER_SEQUENCE"],
    rep_prob     = params["SOUND_REP_PROB"],
    sequence_no  = params["LOC_TRIALS"] * params["LOC_REP"]
    )

# 06. RUN THE EXPERIMENT -----------------------------------------------------------------
control.start(skip_ready_screen=True)

### ------------------ MAIN TASK  ------------------ ##
if main_task_on:
    task_name = "counting"
    
    # Present instructions for the counting task.
    # Wait a minimal time needed to read the instructions.
    instructions_task.present()
    exp.clock.wait(params["WAIT_TIME"])
    keyboard.wait(keys=params['DETECTION_SYMBOL'])

    # Clear instructions and show a blank screen.
    instructions_task.clear_surface()
    blank_canvas.present()

    # Play the soundtrack over the blocks
    for block in range(no_blocks):

        # Correct for zero-indexing
        block_idx = block + 1

        # Select the part of the dataframe relevant for the current block
        df_block = df[df["block_no"] == block_idx]

        # Check if quit key is pressed
        keyboard.check(keys=[misc.constants.K_y])

        # Initialize unique timestamps for logs
        nw = datetime.now()
        ts = int(nw.timestamp())

        # Initialize a dict to track performance on block-level
        trial_performance = {"CORRECT": 0, "INCORRECT": 0}

        # Count sound trials to calculate percentage correct
        no_sound_trials = 0

        # Initialize log for tones
        timDev_log = io.OutputFile(suffix = sesh, directory = f'bids_output_{comboID:003d}')
        timDev_log.write("onset\tduration\ttrial_type\n")
        
        # Initialize log for responses on the frequency deviant counting task
        freqDev_log = io.OutputFile(suffix = sesh, directory = f'bids_output_{comboID:003d}')
        freqDev_log.write("onset\tduration\ttrial_type\tresponse_time\n")

        # Wait for onset of the functional sequence for the main task
        keyboard.wait(keys=[misc.constants.K_s])
        canvas.present()

        # Mark the start of the functional sequence for the main task
        task_start_time = exp.clock.time

        # Wait for 4 's' keys from the scanner to synchronize scanner & script onsets.
        keyboard.wait(keys=[misc.constants.K_s]); keyboard.wait(keys=[misc.constants.K_s])
        keyboard.wait(keys=[misc.constants.K_s]); keyboard.wait(keys=[misc.constants.K_s])

        # Play all tone sequences: trial by trial
        block_start_time = exp.clock.time - task_start_time

        for soundarray, ITI, freq_dev_no, trial_log, time_end in sound_gen.generate_soundtrack(df_block, block_start_time, params["MAX_AMPLITUDE"], params["NUM_HARMONICS"],  params["TONE_DURATION"],  params["HARMONIC_FACTOR"], params["TONE_LOUDNESS"]):

            # Refresh the screen
            canvas.present()

            # Check if quit key is pressed during the trial
            keyboard.check(keys=[misc.constants.K_y])

            # Play the soundtrack
            sd.play(soundarray, samplerate = params["SAMPLE_RATE"])

            # Wait until the end of each trial
            sd.wait()

            # Initialize variables for logging task performance on trial-level
            response, rt = None, None
            max_response_time   = ITI
            response_time_start = clock.time

            # Align the clocks
            time_end_msec = time_end * 1000
            offset_ms   = response_time_start - time_end_msec

            # Check for key presses with max response time equal to ITI
            while clock.time - response_time_start < max_response_time:
                
                # Non-blocking function to check for pressed keys
                keys = keyboard.read_out_buffered_keys()

                # ------ Key-press trials ------ 
                if keys and keys[0] != 115:  # 115 is "s" from scanner sync box
                    response = keys[-1]      # If multiple, we take the last key
                    
                    # Get time of key press
                    key_press_raw = clock.time
                    key_press_corr = key_press_raw - offset_ms
                    key_press_time = key_press_corr / 1000
                    
                    # Get reaction time
                    rt = (clock.time - response_time_start) / 1000

                    # Show feedback
                    perfo_code, trial_performance = give_feedback_freqCount(
                        freq_dev_no, # Actual number of frequency deviants
                        response,    # Participant's count
                        correct,     # Fixation cross for correct responses
                        wrong        # Fixation cross for incorrect responses
                    )

                    # Log with correction of the response for MRI buttons (1 is 0 devs, ...)
                    response = response - 1
                    freqDev_log.write(f"{key_press_time}\t100\t{chr(response)}\t{rt}\n")

            # Update count of sound trials for percentage correct
            if not np.isnan(freq_dev_no):
                no_sound_trials += 1

            # Write relevant info to log for tones
            timDev_log.write(f"{trial_log}")
        
        # Rename the logs according to BIDS standard
        timDev_log.rename(f"sub-{exp.subject:02d}_ses-{sesID:02d}_task-timDev_ts-{ts}_events.tsv")
        freqDev_log.rename(f"sub-{exp.subject:02d}_ses-{sesID:02d}_task-freqDev_ts-{ts}_events.tsv")

        # Save the logs on block level!
        timDev_log.save()
        freqDev_log.save()
        
        # Clearing any key presses
        keyboard.clear()

        # Check if quit key is pressed
        keyboard.check(keys=[misc.constants.K_y])

        # If the task ends before the MRI protocol. Inform the participant to keep calm and remain still.
        blank_canvas.present()
        scanner_text.present()

        # Wait for the MRI QU to end the block by pressing a key unavailable to the participant.
        keyboard.wait(keys=[misc.constants.K_e])

        # When the MRI sequence ends, show text with performance updates
        mainTask_performance = f'Correct: {trial_performance["CORRECT"]}, Wrong: {trial_performance["INCORRECT"]}'
        percent_correct = f'Percentage correct: {round((trial_performance["CORRECT"] / no_sound_trials) * 100)} %'
        mainTask_progress    = f'Blocks completed: {block_idx} / {no_blocks}'
        mainTask_rest = stimuli.TextScreen(
            params["REST_HEADING"], 
            mainTask_performance + f"\n\n{percent_correct}" + f"\n\n{mainTask_progress}" + "\n\n" + params["REST_TEXT"],
            heading_size=params["HEADING_SIZE"],
            heading_colour=params["WHITE"],
            text_size=params["TEXT_SIZE"],
            text_colour=params["WHITE"]
            )
        mainTask_rest.present()

        # Wait for the participant to signal they are rested
        keyboard.wait(keys=params['DETECTION_SYMBOL'])
        mainTask_rest.clear_surface()
        blank_canvas.present()

    # Say thanks at the end of the frequency counting task
    goodbye_task_message = stimuli.TextScreen(
        params["TASK_END_HEADING"] + task_name + " task.",
        params["TASK_END_TEXT"],
        heading_size=params["HEADING_SIZE"],
        heading_colour=params["WHITE"],
        text_size=params["TEXT_SIZE"],
        text_colour=params["WHITE"]
    )
    goodbye_task_message.present()
    exp.clock.wait(params["WAIT_TIME"])
    goodbye_task_message.clear_surface()

### ------------------ LOCALIZER  ------------------ ##
if localizer_on:
    task_name = "repetition detection"
    
    # Decide randomly to start with silence or sound.
    start_with_sound = random.choice([True, False])

    # Present instructions for the repetition detection task.
    # Wait a minimal time needed to read the instructions.
    instructions_loc.present()
    exp.clock.wait(params["WAIT_TIME"])
    keyboard.wait(keys=params['DETECTION_SYMBOL'])

    # Clear and present a blank screen.
    instructions_loc.clear_surface()
    blank_canvas.present()

    # Start the localizer
    loop = 0
    for run in range(params["LOC_REP"]):

        # Initialize the log with unique timestamps
        nw = datetime.now()
        ts = int(nw.timestamp())
        localizer_log = io.OutputFile(suffix = sesh, directory = f'bids_output_{comboID:003d}')
        localizer_log.write("onset\tduration\tstim_file\ttrial_type\tresponse_time\n")
        run_performance = {"H": 0, "M": 0, "CR": 0, "FA": 0}

        # Wait for onset of the MRI sequence. Show the fixation cross at start.
        keyboard.clear()
        keyboard.wait(keys=[misc.constants.K_s])
        canvas.present()

        # Mark the start time of the functional sequence
        run_start_time = exp.clock.time

        # Wait for 4 's' keys from the scanner to synchronize scanner & script onsets
        keyboard.wait(keys=[misc.constants.K_s]); keyboard.wait(keys=[misc.constants.K_s])
        keyboard.wait(keys=[misc.constants.K_s]); keyboard.wait(keys=[misc.constants.K_s])

        # Loop through the trials
        for trial in range(params["LOC_TRIALS"]):

            # Refresh screen
            canvas.present()

            # Check if quit key is pressed
            keyboard.check(keys=[misc.constants.K_y])

            if start_with_sound:

                # Sound part
                sounds_in_sequence = []
                t1, t2 = play_sounds(soundtrack[loop],
                            params["SOUND_DURATION"],
                            exp,
                            canvas,
                            fix_cross,
                            correct,
                            wrong,
                            keyboard,
                            params["DETECTION_SYMBOL"],
                            localizer_log)

                # Refresh the screen
                canvas.present()

                # Silent part
                t3, t4 = play_silence(silence,
                            params["SOUND_DURATION"],
                            exp,
                            params["SOUNDS_PER_SEQUENCE"],
                            keyboard,
                            params["DETECTION_SYMBOL"],
                            localizer_log)
                # print("SOUND FIRST"); compute_durations(params, t1, t2, True); compute_durations(params, t3, t4, True)
                
                # Refresh the screen
                canvas.present()

            else:
                # Silent part
                t5, t6 = play_silence(silence,
                            params["SOUND_DURATION"],
                            exp,
                            params["SOUNDS_PER_SEQUENCE"],
                            keyboard,
                            params["DETECTION_SYMBOL"],
                            localizer_log)

                # Refresh the screen
                canvas.present()

                # Sound part
                sounds_in_sequence = []
                t7, t8 = play_sounds(soundtrack[loop],
                            params["SOUND_DURATION"],
                            exp,
                            canvas,
                            fix_cross,
                            correct,
                            wrong,
                            keyboard,
                            params["DETECTION_SYMBOL"],
                            localizer_log)
                # print("SILENCE FIRST"); compute_durations(params, t5, t6, True); compute_durations(params, t7, t8, True)

                # Refresh the screen
                canvas.present()

            # Update the count
            loop += 1

        # Save the log
        localizer_log.rename(f"sub-{exp.subject:02d}_ses-{sesID:02d}_task-localizer_ts-{ts}_events.tsv")
        localizer_log.save()

        # The experiment ends before the MRI protocol. Inform the participant to keep calm and remain still.
        blank_canvas.present()
        scanner_text.present()

        # Wait for the MRI QU to end the run by pressing a key unavailable to the participant.
        keyboard.wait(keys=[misc.constants.K_e])

        # Give encouragement and performance update on all runs.
        loc_performance = f'Correct: {run_performance["H"]}, Wrong: {run_performance["FA"] + run_performance["M"]}'
        loc_progress    = f'Runs completed: {run+1}/{params["LOC_REP"]}'
        loc_rest = stimuli.TextScreen(
            "PERFORMANCE REPORT", 
            loc_performance + f"\n\n{loc_progress}",
            heading_size=params["HEADING_SIZE"],
            heading_colour=params["WHITE"],
            text_size=params["TEXT_SIZE"],
            text_colour=params["WHITE"]
            )
        loc_rest.present()
        exp.clock.wait(params["WAIT_TIME"])
        loc_rest.clear_surface()
        blank_canvas.present()

    # Say goodbye at the end of the localizer task
    goodbye_task_message = stimuli.TextScreen(
        params["TASK_END_HEADING"] + task_name + " task.",
        params["TASK_END_TEXT"],
        heading_size=params["HEADING_SIZE"],
        heading_colour=params["WHITE"],
        text_size=params["TEXT_SIZE"],
        text_colour=params["WHITE"]
    )
    goodbye_task_message.present()
    exp.clock.wait(params["WAIT_TIME"])
    goodbye_task_message.clear_surface()

# Say thanks & goodbye to the participant
goodbye_exp_message.present()
exp.clock.wait(params["WAIT_TIME"])
goodbye_exp_message.clear_surface()

#7: END EXPERIMENT
control.end()