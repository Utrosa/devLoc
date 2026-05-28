function tapas(homePath, subID, sesID, project, task)
%% Preparation
    % Set up paths to raw DICOM files and BIDS-compliant physio data.
    dcmPath    = sprintf('%s/data_MRI/sourcedata/dicoms/sub-%02d_ses-%02d_%s/', homePath, subID, sesID, project);
    biopacPath = sprintf('%s/data_physio/raw/', homePath);

    % Identify target DICOM images in the DICOM directory.
    % Select folders with 'mm' in the foldername (unique to func)
    localizerDir = dir([dcmPath '*mm*']);
    
    % Exclude non-bold scans: fieldmaps, phase, and full-head scan
    excludeDir   = {'SE', 'Pha', 'FH'};
    localizerDir = localizerDir(~cellfun(@(x) any(contains(x, excludeDir)), {localizerDir.name}));

    % Set up the same func scan labels as in the configuration for raw MRI data
    names = {localizerDir.name};
    n = length(names);

    % Extract the trailing numbers from each name for sorting
    lastTwoNums = zeros(1, n);
    for i = 1:n
        name = names{i};

        % Extract the sequence of digits at the very end of the string
        tokens = regexp(name, '(\d+)$', 'tokens');
        if ~isempty(tokens)
            lastTwoNums(i) = str2double(tokens{1}{1});
        else
            lastTwoNums(i) = i; % Fallback if no number found
        end
    end

    % Get sorting indices (ascending: lowest number first)
    [~, sortIdx] = sort(lastTwoNums);
    
    % Reorder the localizerDir structure based on the sorted indices
    localizerDir = localizerDir(sortIdx);

    % Assign new names sequentially
    for i = 1:n
        if i == n
            % The last element (highest number) gets 'FUNCLOC'
            % Assuming that localizer is the last scan collected
            localizerDir(i).bids_name = 'FUNCLOC';
        else
            % All previous elements get BLOCK1, BLOCK2, etc.
            localizerDir(i).bids_name = sprintf('BLOCK%d', i);
        end
    end
  
    % Remove empty entries (0×0 double, '', missing, etc.)
    bids_names = {localizerDir.bids_name};
    bids_names = bids_names(~cellfun(@isempty, bids_names));

    % Display all the functional scans found for this participant.
    fprintf('%d FUNCTIONAL SEQUENCES FOUND:\n%s\n', numel(bids_names), strjoin(bids_names, ',\n'));

    % Find the first DICOM file per functional scan.
    funcInfo = struct();
    for i = 1:length(bids_names)
        thisVolNum = inf;
        thisDicom = '';

        % Get a list of files in the directory
        tmpDir = dir([localizerDir(i).folder '/' localizerDir(i).name '/']);

        % Loop through all DICOMS per functional scan.
        for dicom = 1:length(tmpDir)
    
            % Skip directories
            if tmpDir(dicom).isdir
                continue;
            end
    
            % Split the DICOM names by the dot of the extension
            ix     = strsplit(tmpDir(dicom).name, '.');
            numVal = str2double(ix{end-1});
    
            % Find the smallest DICOM ID.
            if ~isnan(numVal) && numVal < thisVolNum
                thisVolNum = numVal;
                thisDicom = [tmpDir(dicom).folder '/' tmpDir(dicom).name];
            end
        end
    
        % Check if a DICOM was actually found before proceeding
        if isempty(thisDicom)
            error('No valid DICOM file found in: %s', [localizerDir(i).folder '/' localizerDir(i).name '/']);
        end
    
        % Add the identified first DICOM & other necessary info.
        funcInfo.firstDICOMimage(i) = {thisDicom} ;
    
        % Get info from the first DICOM image
        metaData  = dicominfo(thisDicom);
    
        % The InstanceNumber refers to the number of the file.
        % In the old DICOM format, there is one file / slice.
        % Current file no: InstanceNumber
        % No. of vols/measurements : NumberOfTemporalPositions
        funcInfo.NofVols(i)   = metaData.NumberOfTemporalPositions;
        funcInfo.NofSlices(i) = metaData.NumberOfFrames;

        % No TR is reported in DICOM metadata, so read from filename.
        strTR  = regexp(localizerDir(i).name, 'TR(\d+)', 'tokens');
        funcInfo.TR(i) = str2double(strTR{1}{1});
    end

    clear dicom; clear thisVolNum; clear thisDicom; clear numVal; clear ix;
    clear metaData; clear excludeDir; clear strTR; clear tmpDir; clear i;
    
%% APPLY TAPAS TO THE FUNCTIONAL SEQUENCE/S -------------------------------
    % Specify the functional sequence or loop over a list.
    i = 5; % FUNCLOC acquisition is the last one (5th)
    fprintf('\nAPPLYING TAPAS TO:\n%s', bids_names{i});
    
    %% Structure the BIDS data as TAPAS expects it.
    % Unzip physio file (output from physio2bids)
    gunzip(sprintf( ...
        '%s/data_physio/tmp/sub-%02d_ses-%02d_task-%s_physio_09_125Hz.tsv.gz', ...
        homePath, subID, sesID, task));
    
    % Read the unzipped files into a table.
    origTsvFile = readtable( ...
        sprintf('%s/data_physio/tmp/sub-%02d_ses-%02d_task-%s_physio_09_125Hz.tsv', homePath, subID, sesID, task), ...
        'FileType', 'text', 'Delimiter','\t', 'ReadVariableNames', false);
    
    % Extract the relevant columns for tapas: the column order is known
    % from .json files: time | trigger | puls | resp | CO2 | 02
    phys2bids_data = table2array(origTsvFile);
    raw_trigger = phys2bids_data(:,2);
    cardiac = phys2bids_data(:,3);
    resp = phys2bids_data(:,4);
    
    % Threshold the trigger to contain 0 or 1 values
    % Compute midpoint threshold between low and high trigger levels
    thr = (max(raw_trigger) + min(raw_trigger)) / 2;
    trigger_binary = raw_trigger > thr;
    trigger = [0; diff(trigger_binary) == 1];
    trigger = double(trigger);

    % Arrange the data in the desired order.
    tapas_data = [cardiac resp trigger];
       
    % Overwrite the unzipped data with tapas-compatible file.
    writematrix(tapas_data, ...
    sprintf('%s/data_physio/tmp/sub-%02d_ses-%02d_task-%s_physio_09_125Hz.tsv', homePath, subID, sesID, task), ...
    'FileType','text', ...
    'Delimiter','\t');
    
    % Create the main input structure - PhysIO
    physio = tapas_physio_new();
    
    %% save_dir module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Directory where output files are saved to.
    physio.save_dir = biopacPath;
    
    %% write_BIDS module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    physio.write_bids.bids_step   = 4;
    physio.write_bids.bids_dir    = biopacPath;
    physio.write_bids.bids_prefix = sprintf( ...
        'sub-%02d_ses-%02d_task-%s', subID, sesID, task);
    
    %% log_files module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % General physiological log-file information.
    physio.log_files.vendor = 'BIDS'; % 'Biopac_Mat' or 'Biopac_Txt'

    % Logfiles with cardiac and respiratory data
    % File number 9 corresponds to the functional scan of the localizer for session 2.
    physio.log_files.cardiac =  sprintf('%s/data_physio/tmp/sub-%02d_ses-%02d_task-%s_physio_09_125Hz.tsv', homePath, subID, sesID, task);
    physio.log_files.respiration =  sprintf('%s/data_physio/tmp/sub-%02d_ses-%02d_task-%s_physio_09_125Hz.tsv', homePath, subID, sesID, task);
    
    % Additional file for relative timing information between logfiles 
    % and MRI scans.
    % The time stamp in the DICOM header is on the same time axis as 
    % the time stamp in the physiological log file.
    % If you set this file for "BIDS" input data, it will think the dicom
    % is the json file. Currently, only implemented for Siemens.
    % physio.log_files.scan_timing = funcInfo.firstDICOMimage{i};
    
    % Sampling rate is 125 Hz for all channels except the Trigger
    % channel which is at 250 Hz.
    % phys2bids outputs all channels in either 125 Hz and 250 Hz
    physio.log_files.sampling_interval = 125;

    % For BIDS, Siemens, set to 0.
    physio.log_files.relative_start_acquisition = 0;
    
    % Which scan shall be aligned to which part of the logfile.
    physio.log_files.align_scan  = 'last';
    
    %% scan_timing module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters for sequence timing & synchronization

    % Number of slices per volume
    physio.scan_timing.sqpar.Nslices = funcInfo.NofSlices(i);
    
    % Equals Nslices because we did not trigger with the heart beat.
    physio.scan_timing.sqpar.NslicesPerBeat = funcInfo.NofSlices(i);
    
    % Volume repetition time in seconds
    physio.scan_timing.sqpar.TR = funcInfo.TR(i)/1000;
    
    % Number of full volumes saved: rows in the nifti design matrix.
    physio.scan_timing.sqpar.Nscans = double(funcInfo.NofVols(i));
    
    % Slice whose start determines the adjustment of the regressor
    % timing.
    physio.scan_timing.sqpar.onset_slice = 1;
    
    % Method to determine slice acquisition onset times.
    % 'nominal'           derive slice acquisition timing from sqpar directly
    % 'gradient_log'      derive from logged gradient time courses (Philips)
    % 'scan_timing_log'   uses individual scan timing logfile with time stamps
    %                     specified in log_files.scan_timing
    %                     e.g.,
    %                     *_INFO.log for 'Siemens_Tics' (time stamps for every slice and volume)
    %                     *.dcm (DICOM) for Siemens, is first volume (non-dummy) used in GLM analysis
    %                     *.tsv (3rd column) for BIDS, using the scanner volume trigger onset events
    %                     NOTE: This setting needs a valid filename entered in log_files.scan_timing.
    physio.scan_timing.sync.method = 'scan_timing_log';
    
    %% preproc module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing strategy and parameters for physiological data.
    
    % Measurement modality of input cardiac signal: 'ECG' or 'PPU'
    physio.preproc.cardiac.modality = 'PPU';
    
    % Filter properties for bandpass-filtering of cardiac signal before peak
    % detection, phase extraction, and other physiological traces.
    physio.preproc.cardiac.filter.include = 0; % 1 = YES; 0 = NO
    
    % The initial cardiac pulse selection structure determines how most
    % of the cardiac pulses are detected.
    % 'auto_matched' [default]: auto generation of representative QRS-wave;
    %                           detection via max. auto-correlation with it.
    % 'load_from_logfile': from phys logfile, detected R-peaks of scanner.
    % 'manual': via manually selected QRS-wave for autocorrelations.
    % 'load': from previous manual/auto run.
    physio.preproc.cardiac.initial_cpulse_select.method = 'auto_matched';

    % Maximum allowed physiological heart rate in beats per minute.
    physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
    
    % Peak height threshold in z-scored cardiac waveform to find pulse events.
    physio.preproc.cardiac.posthoc_cpulse_select.min = 0.4;
    
    % The post-hoc cardiac pulse selection structure: If only few (<20) cardiac
    % pulses are missing in a session due to bad signal quality, a manual
    % selection after visual inspection is possible. The results are saved for
    % reproducibility.
    %
    % 'off'     - no manual selection of peaks
    % 'manual'  - pick and save additional peaks manually
    % 'load'    - load previously selected cardiac pulses
    physio.preproc.cardiac.posthoc_cpulse_select.method = 'off';

    % Detecting suspicious outliers
    physio.preproc.cardiac.posthoc_cpulse_select.percentile = 80;
    physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh = 60;
    physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh = 60;

    % [f_min, f_max] frequency interval in Hz of all frequency that should
    % pass the passband filter. Remove high frequency noise and low frequency
    % drifts, but does not distort.
    physio.preproc.respiratory.filter.passband = [0.01, 2.0];
    
    % Whether to remove spikes from the raw respiratory trace using a sliding
    % window median filter.
    physio.preproc.respiratory.despike = false;
    
    %% model module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Derive physiological noise model from preprocessed data. 
    % Models can be combined.
    
    % Unless, we want a session mean (then set to 'all'), no orthogonalisation
    % is needed because our acquisition was NOT triggered to heartbeat.
    physio.model.orthogonalise = 'none';
    
    % True only for RETROICOR model
    physio.model.censor_unreliable_recording_intervals = true;
    
    % Saving the output physio-structure and regressors
    physio.model.output_multiple_regressors = sprintf('sub-%02d_ses-%02d_task-%s_acq-%s_regressors.tsv', subID, sesID, task, string(bids_names{1, i}));
    physio.model.output_physio              = sprintf('sub-%02d_ses-%02d_task-%s_acq-%s_physio.mat', subID, sesID, task, string(bids_names{1, i}));
    
    %%%% RETROICOR Model: Glover et al. 2000. %%%%%%%%%%%%%%%%%%%%%%%%%
    physio.model.retroicor.include = 1;  % 1 = included; 0 = not used
    
    % Natural number, order of cardiac phase Fourier expansion.
    physio.model.retroicor.order.c  = 3;
    
    % Natural number, order of respiratory phase Fourier expansion.
    physio.model.retroicor.order.r  = 4;
    
    % Natural number, order of cardiac-respiratory-phase-interaction 
    % Fourier expansion
    physio.model.retroicor.order.cr = 1;   

    %% verbose module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    physio.verbose.level = -1;

    %% ons_sec module %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    physio.ons_secs.c_scaling = 1;
    physio.ons_secs.r_scaling = 1;
    
    %% Create main regressors by inputing the PhysIO
    [physio_out, R, ons_secs] = tapas_physio_main_create_regressors(physio);
end 