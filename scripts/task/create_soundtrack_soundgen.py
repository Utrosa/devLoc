#! /usr/bin/env python
# Time-stamp: <2026-04-27 m.utrosa@bcbl.eu>
'''
create_soundtrack_soundgen() module generates sounds as defined in csv.
'''
# import warnings # ON PAUSE because it's annoying
import numpy as np
import pandas as pd
from pathlib import Path
import sounddevice as sd
import ast

# TODO: replace set_dbspl() with Jasmin's code for sound normalization
def set_dbspl(sound, dbspl, ref=20e-6):
    """
    Normalize waveform to target dB SPL.

    :param sound : np.array, input waveform
    :param dbspl: float, desired sound level in dB SPL
    :param ref: float, reference pressure (default 20 µPa)
    """

    # Apply dB SPL scaling (RMS based)
    rms = np.sqrt(np.mean(sound**2))
    target_rms = ref * (10 ** (dbspl / 20))
    scale = target_rms / rms
    scaled_sound = sound * scale

    return scaled_sound

class SoundGen:
    def __init__(self, sample_rate, tau):
        """
        Initialize the CreateSound instance.

        :param sample_rate: Sample rate of sounds in Hz.
        :param tau: The ramping window in milliseconds.
        """
        self.sample_rate = sample_rate
        self.tau = tau / 1000 # convert to sec
    
    def sound_maker(self, freq, max_amplitude, num_harmonics, tone_duration, harmonic_factor, dbspl):
        """
        Make a single normalized sound.

        :param freq: Tone frequency in Hz.
        :param max_amplitude: Maximum amplitude to avoid clipping.
        :param num_harmonics: Number of harmonic tones.
        :param tone_duration: Duration of the tone in seconds.
        :param harmonic_factor: Harmonic amplitude decay factor for the tone.
        :param dbspl: Desired dB SPL (loudness) level (cannot change post sound creation).
        
        :return: normalized_sound: an array of audio samples representing a harmonic complex tone.
        """
        
        # Create a time array: each sample represents one event per second 
        t = np.linspace(
            0,                                     # start
            tone_duration,                         # stop
            int(self.sample_rate * tone_duration), # number of samples
            endpoint = False                       # stop is not the last sample
            )
        
        # Initialize a sound array
        sound = np.zeros_like(t)

        # Generate the harmonics
        for k in range(1, num_harmonics + 1):
            harmonic  = np.sin(2 * np.pi * freq * k * t)
            amplitude = max_amplitude * (harmonic_factor ** (k - 1)) / num_harmonics
            sound += amplitude * harmonic

        # Normalize the sound
        normalized_sound = set_dbspl(sound, dbspl)
        
        return normalized_sound

    def sine_ramp(self, sound):
        """ Apply ramping to the start and end of the sound """

        L = int(self.tau * self.sample_rate)
        t = np.linspace(0, L / self.sample_rate, L)
        sine_window = np.sin(np.pi * t / (2 * self.tau)) ** 2  # Sine fade-in

        sound = sound.copy()
        sound[:L] *= sine_window         # Apply fade-in
        sound[-L:] *= sine_window[::-1]  # Apply fade-out

        return sound

    def generate_soundtrack(self, df, current_time, max_amplitude, num_harmonics, tone_duration, harmonic_factor, dbspl):
        """
        Generate tone sequences with timing deviants for current trial.

        :param df: A dataframe with tone sequence parameters with msec as the time unit.
        :param current_time: Current time in the experiment in sec, relative to start of the task.
        :param max_amplitude: Maximum amplitude to avoid clipping.
        :param num_harmonics: Number of harmonic tones.
        :param tone_duration: Duration of the tone in milliseconds.
        :param harmonic_factor: Harmonic amplitude decay factor for the tone.
        :param dbspl: Desired dB SPL (loudness) level (cannot change post sound creation).
        
        :yield: final_sequence: An array of audio samples, representing harmonic a complex tone sequence.
        """
        current_time = current_time / 1000
        current_time = current_time * self.sample_rate

        # Get number of trials per block
        no_trials = len(df["trial_no"].unique())
        
        # Reminder to yourself that we're assuming msec as unit for ISI, ITI, and DEV
        sample_isi = df["isi"].iloc[0] if not df.empty else "N/A"
        sample_iti = df["iti"].iloc[0] if not df.empty else "N/A"
        sample_dev = df["dev"].iloc[0] if not df.empty else "N/A"
        # message = (
        #     "\nAssuming input values are in milliseconds."
        #     "\nThe script converts TONE_DURATION, ISI, ITI, and DEV to seconds. Verify units in the input dataset."
        #     "\nExample of raw (unconverted) values:"
        #     f" TONE_DURATION ({tone_duration}), ISI ({sample_isi}), ITI ({sample_iti}) and DEV ({sample_dev})."
        #     )
        # warnings.warn(message, UserWarning)

        # Convert to sec (only once)
        tone_duration = tone_duration / 1000

        # Loop through all trials in the dataframe
        # Each trial is a linear combination of parameters
        for trial in df.itertuples():

            # Initialize the sequence, log and count of frequency devs.
            sequence = []
            sequence_log = str()
            freq_dev_count = 0

            # Raise error if timing and frequency devs occur on the same tone
            if not np.isnan(trial.dev_loc): # nan for silent trials
                if trial.dev_loc in trial.freq_loc:
                    raise ValueError(
                        f"\nFor trial-{trial.trial_no:02d} block-{trial.block_no:02d}."
                        " Frequency and timing deviations "
                        f"occur on the same tone (idx: {trial.dev_loc})."
                        "\nCheck your input dataframe."
                        " Parameter combinations may be set incorrectly."
                        )

            # Convert isi & iti to sec on every trial
            isi = trial.isi / 1000
            iti = trial.iti / 1000
                        
            # Calculate how many isi & iti events/samples occur per event
            iti_samples  = int(iti * self.sample_rate)
            isi_samples  = int(isi * self.sample_rate)
            tone_samples = int(tone_duration * self.sample_rate)

            # For sound trials: convert dev to sec & get no. of samples
            if not pd.isna(trial.dev_loc):
                dev = trial.dev_abs / 1000
                dev_samples = int(dev * self.sample_rate)

            # Loop through each tone in the sequence
            for i in range(trial.no_tones):

                # Correct for zero indexing
                tone_count = i + 1

                # Initialize for logging
                fDev, fDevLoc = None, None

                # ----------------- Adding TONES ------------------
                ### SILENT trials
                # If the current trial is silent, "dev" is None.
                if pd.isna(trial.dev):
                    sound = np.zeros(tone_samples)

                ### SOUND trials
                # If the current trial has no frequency deviations,
                # the first (and only) element of "freq_dev" is False.
                elif trial.freq_dev[0]:

                    # Generate frequency deviant tone at given locations
                    if trial.freq_dev_no != 0:
                        freq_loc = sorted(trial.freq_loc)
                        if freq_dev_count < trial.freq_dev_no:
                            freq_loc = freq_loc[freq_dev_count]
                            fDevLoc  = freq_loc
                        else:
                            freq_loc = False

                    if tone_count == freq_loc:
                        sound = self.sound_maker(
                            trial.freq_dev[freq_dev_count],
                            max_amplitude,
                            num_harmonics,
                            tone_duration,
                            harmonic_factor,
                            dbspl
                            )
                        fDev = trial.freq_dev[freq_dev_count]
                        
                        # Update frequency dev count
                        freq_dev_count += 1

                    # Generate frequency standard tone at other locations
                    else:
                        sound = self.sound_maker(
                        trial.base_freq,
                        max_amplitude,
                        num_harmonics,
                        tone_duration,
                        harmonic_factor,
                        dbspl
                        )

                # Generate frequency standard tone sequence
                else:
                    sound = self.sound_maker(
                        trial.base_freq,
                        max_amplitude,
                        num_harmonics,
                        tone_duration,
                        harmonic_factor,
                        dbspl
                        )

                # Apply ramp to start and end using the sine_ramp method
                # for sound trials only
                if pd.isna(trial.dev_loc):
                    ramped_sound = np.zeros(tone_samples)
                else:
                    ramped_sound = self.sine_ramp(sound)
                
                # Get tone onset and add to log
                onset_sec = current_time / self.sample_rate

                if not pd.isna(trial.dev):
                    if trial.dev > 0:
                        if not fDev:
                            if tone_count == trial.dev_loc:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_delta-p{int(trial.dev_abs)}ms_tDevLoc-{int(trial.dev_loc)}_type-fStdtDev"
                            else:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_type-fStdtStd"
                        else:
                            if tone_count == trial.dev_loc:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_fDev-{fDev}Hz_fDevLoc-{int(fDevLoc)}_delta-p{int(trial.dev_abs)}ms_tDevLoc-{int(trial.dev_loc)}_type-fDevtDev"
                            else:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_fDev-{fDev}Hz_fDevLoc-{int(fDevLoc)}_type-fDevtStd"
                    elif trial.dev < 0:
                        if not fDev:
                            if tone_count == trial.dev_loc:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_delta-n{int(trial.dev_abs)}ms_tDevLoc-{int(trial.dev_loc)}_type-fStdtDev"
                            else:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_type-fStdtStd"
                        else:
                            if tone_count == trial.dev_loc:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_fDev-{fDev}Hz_fDevLoc-{int(fDevLoc)}_delta-n{int(trial.dev_abs)}ms_tDevLoc-{int(trial.dev_loc)}_type-fDevtDev"
                            else:
                                tone_type = f"fStd-{int(trial.base_freq)}Hz_fDev-{fDev}Hz_fDevLoc-{int(fDevLoc)}_type-fDevtStd"
                    else:
                        if not fDev:
                            tone_type = f"fStd-{int(trial.base_freq)}Hz_type-fStdtStd"
                        else:
                            tone_type = f"fStd-{int(trial.base_freq)}Hz_fDev-{fDev}Hz_fDevLoc-{int(fDevLoc)}_type-fDevtStd"
                else:
                    tone_type = "silence"

                log_format = f"{onset_sec}\t{tone_duration}\t{tone_type}\n"
                sequence_log = sequence_log + log_format
                
                # Add the sound to the sequence
                sequence.append(ramped_sound)
                current_time += tone_samples

                # ----------------- Adding ISI --------------------
                current_isi = isi_samples

                # For sound trials with deviations
                if not pd.isna(trial.dev_type):

                    # Late tones
                    # The ISI before the current tone is longer, ISI after shorter.
                    # "dev_samples" are calculated from absolute timing deviant value
                    if trial.dev_type == "late":
                        if tone_count == (trial.dev_loc - 1):
                            current_isi = isi_samples + dev_samples
                        elif tone_count == trial.dev_loc:
                            current_isi = isi_samples - dev_samples

                    # Early tones
                    # The ISI before the current tone is shorter, ISI after longer.
                    elif trial.dev_type == "early":
                        if tone_count == (trial.dev_loc - 1):
                            current_isi = isi_samples - dev_samples
                        elif tone_count == trial.dev_loc:
                            current_isi = isi_samples + dev_samples
                    
                    # On-time tones
                    else:
                        current_isi = isi_samples

                # Add the ISI
                # Note: there's one less isi in the sequence than tones.
                if tone_count < trial.no_tones:
                    sequence.append(np.zeros(current_isi))
                    current_time += current_isi

            # -------------- Join all segments ----------------
            final_sequence = np.concatenate(sequence)

            # Check that frequency deviants were counted correctly.
            if not pd.isna(trial.freq_dev_no):
                if trial.freq_dev_no != freq_dev_count:
                    raise ValueError(
                        f"Counted more/less frequency deviants ({freq_dev_count}) "
                        f"than specified ({trial.freq_dev_no}) "
                        f"for trial {trial.trial_no} in block {trial.block_no}.")

            # Yield the tone sequence, ITI (an array of zeros),
            # the number of frequency deviants, and the sequence log.
            end_time = current_time / self.sample_rate
            yield final_sequence, trial.iti, trial.freq_dev_no, sequence_log, end_time

            # Clear the list for the next iteration (memory <3)
            sequence.clear()

            # Add the iti to current time
            current_time += iti_samples

# TEST: example usage -----------------------------------------------------------------------------
if __name__ == "__main__":
    
    # Set the parameters
    sesID = 27
    params = {
        "PROJECT_ROOT"    : "/home/mutrosa/Documents/projects/auditory_paradigms/detection_accuracy/test",  
        "TONE_LOUDNESS"   : 75,     # dB SPL
        "TONE_DURATION"   : 50,     # msec
        "NUM_HARMONICS"   : 10,     # Number of harmonics
        "HARMONIC_FACTOR" : 0.8,    # Harmonic amplitude decay factor
        "MAX_AMPLITUDE"   : 1.14,   # Defined through a simulation
        "SAMPLE_RATE"     : 48000,  # Hz
        "TAU"             : 5,      # Ramping window in msec
        }

    # Load the trial parameters from csv
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

    # Initialize the class
    sound_gen = SoundGen(params["SAMPLE_RATE"], params["TAU"])

    # Play the soundtrack over the blocks
    for i in range(no_blocks):

        # Correct for zero-indexing
        block_idx = i + 1

        # Select only the part of the dataframe relevant for the current trial
        df_block = df[df["block_no"] == block_idx]

        # Generate the soundtrack of the experimental session
        for soundtrack, iti, freq, log in sound_gen.generate_soundtrack(
           df_block,
           block_idx, # This should be the block's start time
           params["MAX_AMPLITUDE"],
           params["NUM_HARMONICS"], 
           params["TONE_DURATION"], 
           params["HARMONIC_FACTOR"],
           params["TONE_LOUDNESS"]
           ):
            sd.play(soundtrack, samplerate = params["SAMPLE_RATE"])
            sd.wait()
            # Here you have to add code to wait for the ITI
            # or output ITI samples and play that as silence with sounddevice