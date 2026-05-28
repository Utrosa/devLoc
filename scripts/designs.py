def localizer(logfilepath):
    """
    Parse logfiles into design matrix in NiPype Bunch format.

    Parameters:
        logfilepaths (list): List of file paths to logfiles.

    Returns:
        list: A list of Bunch objects containing design information.
    """
    import csv
    from nipype.interfaces.base import Bunch
    
    # Get info on stimuli onset, duration and key presses.
    sounds, silences, keypress, sound_prev = [], [], [], []
    with open(logfilepath, 'r') as logfile:
        
        next(logfile) # Skip header row
        
        # Determine the delimiter of the logfiles automatically
        # BIDS-standard assumes tab-separated .event files ;)
        sample  = logfile.read(3000); logfile.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters = [";", "\t" , ","])
        logTsv  = csv.reader(logfile, delimiter="\t")
        
        next(logTsv) # Skip header row again
        for line in logTsv:
            event     = {'onset': float(line[0]), 'duration': float(line[1])}
            stim_file = line[2]

            # Silences
            if stim_file == 'null_event.wav':
                silences.append(event)

            # Sounds with key press during
            elif stim_file.startswith('s3'):
                if stim_file != sound_prev:
                    sounds.append(event)
                    sound_prev = stim_file
                else:
                    if line[4] != 'n/a':
                        keypress.append(event)

            # Sounds with key press after
            elif stim_file == 'n/a':
                if line[4] != 'n/a':
                    keypress.append(event)
            else:
                print('WARNING: Skipping unrecognised line "{}"'.format(line))

    # Incorporate into design info
    conditions = ['sound', 'silence', 'keypress']
    onsets     = [[on['onset'] for on in cond] for cond in [sounds, silences, keypress]]
    durations  = [[du['duration'] for du in cond] for cond in [sounds, silences, keypress]]
    design_info = Bunch(conditions = conditions,
                        onsets     = onsets,
                        durations  = durations)
    return design_info

def timDev(logfilepaths, pooling):
    """
    Parse logfiles into design matrix for the 'timDev' paradigm. 
    Timing deviancy conditions can be taken as absolute or relative values.
    Zero is not included as a timing deviancy condition.

    Parameters:
        logfilepath (str): Path to the log file.
        pooling: If True, timing deviants are pooled as abolute values. 
                 If False, separate conditions for negative and positive values.

    Returns:
        list: A list of Bunch objects containing conditions, onsets, and durations.
    """
    import csv
    from nipype.interfaces.base import Bunch

    design_info_list = []

    for logfilepath in logfilepaths:
        
        # Get info on stimuli onset, duration and key presses from events.tsv file
        # Initialize a dictionary to store that info per timing deviation
        events_by_dev = {}

        with open(logfilepath, 'r') as logfile:
            
            next(logfile)  # Skip header row

            # Auto-detect delimiter (should be tab)
            sample  = logfile.read(3000); logfile.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=[";", "\t", ","])
            
            # Read the logfile
            logTsv  = csv.reader(logfile, dialect)
            next(logTsv)  # Skip header again

            for line in logTsv:

                # Get events
                onset = float(line[0])
                duration = float(line[1])
                event = {'onset': onset, 'duration': duration}
                
                # Get stimulus type
                deviation_str = line[2]
                
                # Initialize deviation to a default value (e.g., None)
                deviation = None

                # Does the current row correspond to a time deviant tone?
                if "delta" in deviation_str:

                    # Strip to get the delta
                    delta_str = deviation_str.split("delta-")[1]
                    delta = delta_str.split("ms")[0]
                    
                    # Figure out the direction: positive or negative delta?
                    if "p" in delta:
                        pD = delta.strip("p")
                        deviation = float(pD)
                    elif "n" in delta:
                        nD = delta.strip("n")
                        if pooling:
                            deviation = float(nD)
                        else:
                            deviation = -float(nD)
                    else:
                        deviation = None

                # Initialize list for this deviation
                if deviation is not None:
                    if deviation not in events_by_dev:
                        events_by_dev[deviation] = []
                
                    events_by_dev[deviation].append(event)

        # Sort deviations from negative to positive
        # Important to ensure consistent order in the conditions list
        sorted_deviations = sorted(events_by_dev.keys())

        # Create conditions (a list of strings)
        conditions = [str(i) for i in sorted_deviations]

        # Extract onsets and durations in the same order as conditions
        onsets = []
        durations = []
        for dev in sorted_deviations:
            onsets.append([e['onset'] for e in events_by_dev[dev]])
            durations.append([e['duration'] for e in events_by_dev[dev]])

        design_info = Bunch(
            conditions=conditions,
            onsets=onsets,
            durations=durations
        )

        # Append to list
        design_info_list.append(design_info)
    
    return design_info_list
