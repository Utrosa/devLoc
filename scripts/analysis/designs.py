#! /usr/bin/env python
# Time-stamp: <10-09-2026 m.utrosa@bcbl.eu>
'''
Defines functions to create SPM design matrices,
needed for NiPype workflow. These matrices have
a format of a Bunch object.

Zero is not included as a timing deviancy condition!

The difference between timDev and timDevCat:
- timDev either takes absolute or nominal values of timing
  deviations, as defined in the experimental task
- timDevCat allows user-specified categorization of timing
  deviations, so making theory-informed grouping (less precise)

TODO: systemize arguments -- can they all have the same name, so it's easier
to adapt the analysis for a different design?
TODO: why do some designs and (timDevCat, localizer) take only one log file
as input and other a list of log paths?
TODO: extend regressors for freqDev (not for my projects)

'''
def localizer(logfilepath):
    """
    Parse a logfile into design matrix in Nipype Bunch format.

    Parameters:
        logfilepath (str): File path to a logfile.

    Returns:
        list: A list of Bunch objects containing design information.
    """
    import csv, warnings
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
    
    # Print an example of the Bunch conditions to terminal
    warnings.warn(
        "\nThe conditions in the bunch of the design `localizer` are:"
        f"\n{design_info.conditions}\n"
    )

    return design_info

def timDev(logfilepaths, absolute):
    """
    Parse logfiles into design matrix for the 'timDev' paradigm. 
    Timing deviancy conditions can be taken as absolute or relative values.
    Zero is not included as a timing deviancy condition.

    Parameters:
        logfilepath (str): Path to the log file.
        absolute: If True, timing deviants are pooled as absolute values. 
                  If False, separate conditions for negative and positive values.

    Returns:
        list: A list of Bunch objects containing conditions, onsets, and durations.
              The list is flattened: no nesting per runs, session, subjects ...
    """
    import csv, warnings
    from nipype.interfaces.base import Bunch

    # Loop through the logs
    design_info_list = []
    for logfilepath in logfilepaths:
        
        # Get info on stimuli onset, duration and trial type from events.tsv file
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
                        deviation = int(pD)
                    elif "n" in delta:
                        nD = delta.strip("n")
                        if absolute:
                            deviation = int(nD)
                        else:
                            deviation = -int(nD)
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

        # Create the bunch object
        design_info = Bunch(
            conditions=conditions,
            onsets=onsets,
            durations=durations
        )

        # Append to list
        design_info_list.append(design_info)
    

    # Print an example of the Bunch conditions to terminal
    warnings.warn(
        "\nThe conditions in the bunch of the design `timDev` are:"
        f"\n{design_info_list[0].conditions}\n"
    )

    return design_info_list

def timDevCat(time_logs, binary, absolute, groups):
    """
    Parse logfiles into design matrix for the 'timDev' paradigm. 
    Timing deviancy conditions can be taken as absolute or relative values.
    Zero is not included as a timing deviancy condition.

    Arguments:
        time_logs (list):
            Paths to the log file.
        binary (bool):
            If True, all timing deviants are grouped into a single condition,
            ignoring magnitude and direction.
        absolute (bool):
            If True, timing deviants are pooled as abolute values. 
            If False, separate conditions for negative and positive values.
        groups (dict):
            A sorted dictionary of upper bounds (keys) and names for the timing 
            deviants groups they create (values). Defaults to False (no grouping).

    Returns:
        list: A Bunch object containing conditions, onsets, and durations.
    """
    import csv, warnings
    from nipype.interfaces.base import Bunch
    from utils import find_dev_group
    
    # Get group values and names, if grouping applies to timing deviants 
    if groups:
        group_values = list(groups.keys())
        group_names  = list(groups.values())

    # Loop through the logs
    design_info_list = []
    for time_log in time_logs:

        # Initialize a dictionary to store info per timing deviation
        events_by_dev = {}
        
        # Read info on stimuli onset, duration and trial type from events.tsv file
        with open(time_log, 'r') as logfile:

            # Skip header row
            next(logfile)

            # Auto-detect delimiter (should be tab)
            sample = logfile.read(3000)
            logfile.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=[";", "\t", ","])
            
            # Read the logfile
            logTsv  = csv.reader(logfile, dialect)
            next(logTsv)  # Skip header again

            for line in logTsv:

                # Get event's onset and duration
                onset = float(line[0])
                duration = float(line[1])
                event = {'onset': onset, 'duration': duration}
                
                # Get stimulus type
                deviation_str = line[2]
                
                # Initialize deviation to a default value (e.g., None)
                deviation = None

                # Does the current row correspond to a time deviant tone?
                if "delta" in deviation_str:
                    if binary:
                        deviation = "timDev"
                    else:
                        # Strip to get the delta
                        delta_str = deviation_str.split("delta-")[1]
                        delta = delta_str.split("ms")[0]
                        
                        # Figure out the direction: positive or negative delta?
                        # Positive deltas
                        if "p" in delta:
                            pD = delta.strip("p")

                            if groups:
                                deviation = find_dev_group(int(pD), groups)
                            else:
                                deviation = int(pD)
                        
                        # Negative deltas                   
                        elif "n" in delta:
                            nD = delta.strip("n")

                            if absolute:
                                if groups:
                                    deviation = find_dev_group(int(nD), groups)
                                else:
                                    deviation = int(nD)
                            else:
                                if groups:
                                    deviation = find_dev_group(-int(nD), groups)
                                else:
                                    deviation = -int(nD)
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

        # Create the bunch object
        design_info = Bunch(
            conditions=conditions,
            onsets=onsets,
            durations=durations
        )

        # Append to list
        design_info_list.append(design_info)

    # Print an example of the Bunch conditions to terminal
    warnings.warn(
        "\nThe conditions in the bunch of the design `timDevCat` are:"
        f"\n{design_info_list[0].conditions}\n"
    )

    return design_info

def freqDev(logfilepaths):
    """
    Parse logfiles into design matrix for the 'freqDev' paradigm. 
    Frequency conditions are either 0 (standard) or 1 (deviant).

    Parameters:
        logfilepath (str): Path to the log file with task == "freqDev"

    Returns:
        list: A list of Bunch objects containing conditions, onsets, and durations.
              The list is flattened: no nesting per runs, session, subjects ...
    """
    import csv, warnings
    from nipype.interfaces.base import Bunch

    design_info_list = []

    for logfilepath in logfilepaths:
        
        # Get info on stimuli onset, duration and trial type from events.tsv file
        # Initialize a dictionary to store that info per timing deviation
        events_by_freq = {'0': [], '1': []}

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
                
                # Get trial type
                frequency_str = line[2]

                # Is it a frequency deviant or not?
                if frequency_str == "0":
                    frequency = "0"
                else:
                    frequency = "1"

                # Initialize a list for this frequency and append
                events_by_freq[frequency].append(event)

        # Create conditions
        conditions = [i for i in events_by_freq.keys()]

        # Extract onsets and durations in the same order as conditions
        onsets = []
        durations = []
        for freq in events_by_freq.keys():
            onsets.append([e['onset'] for e in events_by_freq[freq]])
            durations.append([e['duration'] for e in events_by_freq[freq]])

        design_info = Bunch(
            conditions=conditions,
            onsets=onsets,
            durations=durations
        )

        # Append to list
        design_info_list.append(design_info)
        
    # Print an example of the Bunch conditions to terminal
    warnings.warn(
        "\nThe conditions in the bunch of the design `freqDev` are:"
        f"\n{design_info_list[0].conditions}\n"
    )

    return design_info_list

def timfreqDev(time_log, time_binary, time_abs, time_groups, add_freqDev):
    """
    Joins the events of timDev and freqDev tasks from timDev log file. 

    Parameters:
        time_log:    path to the log files for timDev task.
        time_binary: If True, all timing deviants are grouped into a single condition,
                     ignoring magnitude and direction. Defaults to False.
        time_abs:    If True, timing deviants are pooled as abolute values. 
                     If False, separate conditions for negative and positive values.
        time_groups: dict/bool, a sorted dictionary of upper bounds (keys) and names for the timing 
                     deviants groups they create (values). Default to False (no grouping).
        add_freqDev: If True (str), adding freqDev Bunch. If False, returning only timDev Bunch.
    Returns:
        list: A Bunch object containing conditions, onsets, and durations.
    """
    import csv, warnings
    from nipype.interfaces.base import Bunch
    
    # Nest timDevCat function which takes a single log file as input
    def timDevCat1(time_log, binary, absolute, groups):
        """
        Parse logfiles into design matrix for the 'timDev' paradigm. 
        Timing deviancy conditions can be taken as absolute or relative values.
        Zero is not included as a timing deviancy condition.

        Parameters:
            time_log (str): Path to the log file.
            binary:    If True, all timing deviants are grouped into a single condition,
                       ignoring magnitude and direction. 
            absolute:  If True, timing deviants are pooled as abolute values. 
                       If False, separate conditions for negative and positive values.
            groups:    A sorted dictionary of upper bounds (keys) and names for the timing 
                       deviants groups they create (values). False means no grouping.
        Returns:
            list: A Bunch object containing conditions, onsets, and durations.
        """
        from utils import find_dev_group
        # Initialize a dictionary to store that info per timing deviation
        events_by_dev = {}

        # Get group values and names, if grouping applies to timing deviants 
        if groups:
            group_values = list(groups.keys())
            group_names  = list(groups.values())

        # Read info on stimuli onset, duration and trial type from events.tsv file
        with open(time_log, 'r') as logfile:

            # Skip header row
            next(logfile)

            # Auto-detect delimiter (should be tab)
            sample = logfile.read(3000)
            logfile.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=[";", "\t", ","])
            
            # Read the logfile
            logTsv  = csv.reader(logfile, dialect)
            next(logTsv)  # Skip header again

            for line in logTsv:

                # Get event's onset and duration
                onset = float(line[0])
                duration = float(line[1])
                event = {'onset': onset, 'duration': duration}
                
                # Get stimulus type
                deviation_str = line[2]
                
                # Initialize deviation to a default value (e.g., None)
                deviation = None

                # Does the current row correspond to a time deviant tone?
                if "delta" in deviation_str:
                    if binary:
                        deviation = "timDev"
                    else:
                        # Strip to get the delta
                        delta_str = deviation_str.split("delta-")[1]
                        delta = delta_str.split("ms")[0]
                        
                        # Figure out the direction: positive or negative delta?
                        # Positive deltas
                        if "p" in delta:
                            pD = delta.strip("p")

                            if groups:
                                deviation = find_dev_group(int(pD), groups)
                            else:
                                deviation = int(pD)
                        
                        # Negative deltas                   
                        elif "n" in delta:
                            nD = delta.strip("n")

                            if absolute:
                                if groups:
                                    deviation = find_dev_group(int(nD), groups)
                                else:
                                    deviation = int(nD)
                            else:
                                if groups:
                                    deviation = find_dev_group(-int(nD), groups)
                                else:
                                    deviation = -int(nD)
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
        return design_info
    
    # Create timing deviancy Bunch 
    time_bunch = timDevCat1(time_log, time_groups, time_binary, time_abs)

    # Create frequency deviancy Nunch
    if add_freqDev:
        
        # Get info on stimuli onset, duration and trial type from events.tsv file
        # Initialize a dictionary to store that info per timing deviation
        events_by_freq = {'freqDev': []}

        with open(time_log, 'r') as logfile:
            
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
                
                # Get trial type
                frequency_str = line[2]

                # Is it a frequency deviant or not?
                if "type-fDev" in frequency_str:
                    frequency = "freqDev"
                else:
                    frequency = False

                # Add only frequency deviants as events
                if frequency:
                    events_by_freq[frequency].append(event)

        # Create conditions
        conditions = [i for i in events_by_freq.keys()]

        # Extract onsets and durations in the same order as conditions
        onsets = []
        durations = []
        for freq in events_by_freq.keys():
            onsets.append([e['onset'] for e in events_by_freq[freq]])
            durations.append([e['duration'] for e in events_by_freq[freq]])

        freq_bunch = Bunch(
            conditions=conditions,
            onsets=onsets,
            durations=durations
        )

        # Join the timing and frequency deviancy Bunch objects
        timfreq_bunch = [Bunch(
            conditions=time_bunch.conditions + freq_bunch.conditions,
            onsets=time_bunch.onsets + freq_bunch.onsets,
            durations=time_bunch.durations + freq_bunch.durations
        )]
    else:
        timfreq_bunch = time_bunch

    # Print an example of the Bunch conditions to terminal
    warnings.warn(
        "\nThe conditions in the bunch of the design `timfreqDev` are:"
        f"\n{timfreq_bunch[0].conditions}\n"
    )

    return timfreq_bunch