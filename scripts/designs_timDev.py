def timDev(logfilepaths):
    """
    Parse logfiles into design matrix for the 'timDev' paradigm with 21 conditions.

    Parameters:
        logfilepath (str): Path to the log file.

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
