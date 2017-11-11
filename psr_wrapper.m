function psr_wrapper(parameters) 

% PSR_WRAPPER - Wrapper for processing data with PASER.
% This function loads raw extracellular data, performs all data processing
% steps and saves the results to a MAT file. The processing pipeline
% consists of spike, local field potential and stimulus detection, as well
% as spike sorting.
%
% Syntax:  psr_wrapper(parameters)
%
% Inputs:
%    parameters - See PSR_PARAMETER_DEFAULT and PSR_BATCH_PROCESSING
%
% Outputs:
%    One or more MAT files. See toolbox README for further details.
%
% See also: PSR_BATCH_PROCESSING

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

% fig = figure; set(gcf,'position',get(0,'screensize')); % TEMP

%% Check input

if (~isfield(parameters,'subject')); subject = [];
else,                                subject = parameters.subject;
end

if (~isfield(parameters,'session')); session = [];
else,                                session = parameters.session;
end

if (~isfield(parameters,'loadPathSub')); loadPath = [];
else,                                    loadPath = parameters.loadPathSub;
end

if (~isfield(parameters,'savePathSub')); savePath = [];
else,                                    savePath = parameters.savePathSub;
end

%% Check if paths have correct format

nTrials = length(loadPath);
for iTrial = 1:nTrials
    if (loadPath{iTrial}(end) ~= '\' && ~isempty(loadPath{iTrial}))
        loadPath{iTrial} = [loadPath{iTrial}, '\'];
    end
end

if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

%% Convert initial parameters

temp.path    = parameters.path;
parameters   = rmfield(parameters,'path');
temp.general = parameters; 
parameters   = temp; clear temp;

%% Load processing parameters

if (~isfield(parameters,'configPath')); psr_parameter_default; % Load parameters 
else,                                   run(parameters.configPath);
end
parameters = orderfields(parameters);

%% Constants
nElectrodes  = parameters.general.nelectrodes;
pattern      = parameters.general.filepattern;
ext          = parameters.general.extension;
method       = lower(parameters.sorting.method);
precision    = 10^round(parameters.general.precision);
tempFileName = 'Temp_vars';

%% Check what parts of pipeline need to be run

tf = strcmp(method,{'fmm','ums','ops'});
tf = max(tf(:));
tf = tf & parameters.process.spikes;
if (tf); SPIKE_DETECTION = 1; % Spike detection through thresholding
else;    SPIKE_DETECTION = 0; % Use template matching technique
end

tf = strcmp(method,{'kst','cbp','ops'});
tf = max(tf(:)) |  parameters.spikes.artifacts_removal;
if (tf); SAVE_RAW_DATA = 1; % Save temporary raw data
else;    SAVE_RAW_DATA = 0;
end

%% Check if TEMP files exist in save folder

filesTemp  = cell(0,0);
files  = dir([savePath 'Temp_Probe_*_Trial_*.mat']);
files  = char(files.name);
nfiles = size(files,1);
for iFile = 1:nfiles
    filename = [savePath files(iFile,:)];
    probe = extractStringFromPath(filename,'Probe_(\d*)_');
    trial = extractStringFromPath(filename,'Trial_(\d*).');
    filesTemp{probe,trial} = filename;
end
fileTemp = dir([savePath tempFileName '.mat']);
fileTemp = char(fileTemp.name);

%% Start data processing

if (isempty(filesTemp) || isempty(fileTemp))
    
    %% Extract stimulus conditions from folder name
    
    stimuliConditions = zeros(nTrials,1);
    
    for iTrial = 1:nTrials
        filepath = loadPath{iTrial};
        stimuliConditions(iTrial) = extractStringFromPath(filepath,['_(\d*)' parameters.general.trialpattern]);
    end
    
    %% Find raw data files
    
    filesUnsorted = cell(nTrials,1);
    for iTrial = 1:nTrials
        files = dir([loadPath{iTrial} '\*' pattern '*' ext]);
        if (size(files,1) == 0); return; end
        filesUnsorted{iTrial} = char(files.name);
    end
    
    %% Sort raw data files in correct chronological order 
    
    filesRawAll = cell(nTrials,1);
    
    for iTrial = 1:nTrials
        nFiles = length(filesUnsorted{iTrial}(:,1));
        filesRaw = cell(nFiles,1);
        for iFile = 1:nFiles
            filename   = filesUnsorted{iTrial}(iFile,:);
            filename   = strtrim(filename);
            k          = strfind(filename,pattern) + length(pattern);
            [~,name,~] = fileparts(filename); % remove extension
            name       = name(k:end); % take string after pattern
            k          = strfind(name,'_'); % check if underscore is present
            if (~isempty(k)); name = name(1:k-1); end % take number between pattern and underscore
            id         = str2double(name); % convert to array index
            filesRaw{id} = filename;
        end
        filesRaw = filesRaw(~cellfun('isempty',filesRaw)); % remove empty cells
        filesRawAll{iTrial} = filesRaw;
    end
    
    %% STIMULUS ONSET DETECTION 
    
    % Different detection methods for stimulus onset times. This is heavily
    % dependent on the experimental paradigm.
        
    %%%% ACTIVE VS. PASSIVE
    
    if (parameters.process.active_passive)
        
        stimTimes = cell(nTrials,2); % stimulus onset times
        nSessions = size(parameters.general.session,2);
        sessions  = cell(nSessions,1);
        for iSession = 1:nSessions % Convert to lower case
            sessions{iSession} = lower(parameters.general.session{iSession});
        end
        
        % Find active and passive trials
        k = find(~cellfun(@isempty,strfind(sessions,'passive')));
        k = (parameters.general.sessionTrial == k);
        loadPathActive  = loadPath(~k);
        loadPathPassive = loadPath( k);
        
        % Camera onset times
    
        % For active trials
        stimTimes(~k,1) = psr_cam_detection(loadPathActive);        
    
        % Magnetic field artifacts [MFA]
    
        % For passive trials
        parameters.general.stims = stimuliConditions(k);
        stimTimes(k,:) = psr_mfa_detection(loadPathPassive,parameters);
    
    end
    
    %% SPIKE + LFP DETECTION
        
    nProbes = zeros(nTrials,1);
    parameters.Fr = parameters.lfp.rsfactor * parameters.lfp.bp_upper;
    
    for iTrial = 1:nTrials
        
        % Check for existing files
        
        if (iTrial <= size(filesTemp,2))
            filesRaw = cell2mat(filesTemp(:,iTrial));
            nFiles   = size(filesRaw,1);
            nProbes  = size(filesRawAll{iTrial},1) / parameters.general.nelectrodes;
            if (nFiles == nProbes); continue; end
        end
        
        % Initialize
        
        filesRaw        = filesRawAll{iTrial};
        nFiles          = length(filesRaw);
        nProbes(iTrial) = nFiles / parameters.general.nelectrodes;
        stimTimesTrial  = stimTimes(iTrial,:);
        stim            = stimuliConditions(iTrial);
        
        for iProbe = 1:nProbes(iTrial)
            
            for iElectrode = 1:nElectrodes
                
                % Load Open-Ephys data
                
                iFile = (iProbe - 1) * nElectrodes + iElectrode;
                
                filename = [loadPath{iTrial} filesRaw{iFile}];
                filename = strtrim(filename);
                
                try % Load CONTINUOUS files [microvolts]
                    disp(['Loading ' filename '...']);
                    [data_channel_raw, timestamps, info] = load_open_ephys_data_faster(filename);
                    timestamps = timestamps / info.header.sampleRate;
                catch
                    [data_channel_raw, timestamps, info] = load_open_ephys_data(filename);
                end
                
                % Store variables
                if (iElectrode == 1) % Initialize
                    Fs_array      = zeros(nElectrodes,1);
                    sLength_array = zeros(nElectrodes,1);
                end
                
                Fs = info.header.sampleRate;
                Fs_array(iElectrode) = Fs; % Sampling frequencies should be equal
                sLength_array(iElectrode) = length(data_channel_raw); 
                
                data_channel_raw = psr_artifact_fft(data_channel_raw,parameters,Fs);
                
                if (SPIKE_DETECTION || SAVE_RAW_DATA)
                    
                    %% Initialize arrays
                    if (iElectrode == 1)
                        data_probe_sst = zeros(nElectrodes,length(data_channel_raw),'single');
                    end
                    
                    %% BAND-PASS FILTERING SPIKE SORTING
                    
                    [B,A]        = butter(parameters.spikes.bp_order,[parameters.spikes.bp_lower parameters.spikes.bp_upper]/(Fs/2),'bandpass');
                    data_channel = filtfilt(B,A,data_channel_raw); % Zero-phase digital filtering
                    data_channel = single(data_channel); % convert to single type
                    
                    data_probe_sst(iElectrode,:) = detrend(data_channel, 'linear');
                    
                end
                
                %% LOCAL FIELD POTENTIAL
                if (parameters.process.lfp)
                    % Resample raw signal

                    timestamps        = timestamps - timestamps(1);
                    [data,timestamps] = resample(data_channel_raw,timestamps,parameters.Fr); % Initial resample to avoid errors in ft_preprocessing

                    if (iElectrode == 1)
                        data_probe_lfp = zeros(nElectrodes,length(data),'single');
                        time_probe_lfp = zeros(nElectrodes,length(data),'single');
                    end

                    time_probe_lfp(iElectrode,:) = timestamps;
                    data_probe_lfp(iElectrode,:) = data;
                end
            end
            
            clear data_channel data_channel_raw
            
            % Check data
            
            flag = sum(Fs_array ~= Fs_array(1)) ~= 0;
            if (flag); disp('Sampling frequency mismatch'); continue;
            else; Fs = Fs_array(1);
            end
            
            flag = sum(sLength_array ~= sLength_array(1)) ~= 0;
            if (flag); disp('Data length mismatch'); continue;
            else; sLength = sLength_array(1);
            end
            
            % SPIKE DETECTION
            
            spikes = [];
            spikes.Fs = Fs;  % Hz, sampling rate of spike data
            
            if (SPIKE_DETECTION)
                if (stim > 0 && parameters.mfa_combine) % magnetic field artifact detection
                    spikes = ss_mfa_detection_raw(data_probe_sst,spikes,iStart);
                end

                spikes = psr_sst_detect(data_probe_sst',spikes);
            end
            
            if (parameters.process.spikes)
                if (SAVE_RAW_DATA)
                    spikes.data = int16(precision * data_probe_sst);
                end
            end
            
            % LOCAL FIELD POTENTIAL
            
            freq = [];
            if (parameters.process.lfp && ~isempty(stimTimesTrial{1})) % TODO: REMOVE STIMULUS ONSET TIME CONDITION 
                timestamps = mean(time_probe_lfp,1);
                [data,artifacts] = psr_lfp_artifact_removal(data_probe_lfp,parameters);
                [freq,parameters] = psr_lfp_wrapper(data,timestamps,stimTimesTrial{1},parameters);
                freq.artifacts = artifacts;
            end
            
            % Save
            
            % Create metadata variable: to be expanded in future versions
            % with data from the electronics notebook
    
            metadata = [];
            metadata.subject   = subject;
            metadata.session   = session;
            metadata.dir       = loadPath;
            metadata.stimtimes = stimTimesTrial;
            metadata.duration  = (sLength - 1) / Fs;
            metadata.stimulus  = stim;
            metadata.probe     = iProbe;
            
            if (parameters.process.spikes)
                % Save temporary MAT file
                filename = [savePath 'Temp_Probe_' num2str(iProbe,'%02d') '_Trial_' num2str(iTrial,'%02d') '.mat'];
                filesTemp{iProbe,iTrial} = filename;
                save(filename,'spikes','metadata','freq','parameters');
            else
                % Save LFP output
                saveFile(spikes,[],metadata,freq,parameters,savePath);
            end
        end
        
        %% Filter spiking data across all probes
                        
        if (parameters.process.spikes)
            
            filesTrial = filesTemp(:,iTrial);
            
            if (stim > 0 && parameters.mfa_combine)
                stimTimesTrial = ss_mfa_combine(filesTrial,stimTimesTrial); % Magnetic field artifact combination
                psr_mfa_save(filesTrial,stimTimesTrial)
            end
        end
    end
    
    if (parameters.process.spikes)
        keepvars = {...
            'parameters',        ...
            'filesTemp',         ...
            'nspikes',           ...
            'stimuliConditions', ...
            'metadata',          ...
            'loadPath',          ...
            'savePath',          ...
            'tempFileName'};
        clearvars('-except', keepvars{:});
        save([savePath tempFileName]); % Temporarily save all workspace variables. Useful when continuing after early stopping
    end
else
    load([savePath tempFileName])
end

%% Spike sorting

if (parameters.process.spikes)
    
    
    method     = lower(parameters.sorting.method);
    nProbes    = size(filesTemp,1);
    nTrials    = size(filesTemp,2);
    
    if (~exist('filesSaved','var')) % Initialize array if not loaded
        filesSaved = cell(nProbes,1); % Output files 
    end
    
    for iProbe = 1:nProbes % Do clustering per probe across all stimulus conditions
        
        % Check if output file + temp file (after sorting) exist
        
        if (exist(filesSaved{iProbe},'file')); continue; end 
        
        % Check for missing temporary files
        
        MISSING_FILE = false;
        for iTrial = 1:nTrials
            if (~exist(filesTemp{iProbe,iTrial},'file')); MISSING_FILE = true; break; end
        end
        if (MISSING_FILE) % If not all files are present
            for iTrial = 1:nTrials
                filename = filesTemp{iProbe,iTrial};
                if (exist(filename,'file')); delete(filename); end % Delete temporary spikes file
            end
            disp(['Skipping probe ' num2str(iProbe) '...']);
            continue; % Move on to next probe
        end
        
        disp(['Spike sorting probe ' num2str(iProbe) '...']);
        
        % Initialize
        
        switch method
            case 'cbp'
                dataProbe      = []; %
                dataProbe.data = [];
            case {'fmm','ops','ums'}
                spikesAll = []; % structure to contain all 'spikes' structures from all trials
            case 'kst'
                dataProbe = [];
        end
        
        %% MERGE DATA ACROSS TRIALS
        
        for iTrial = 1:nTrials
            
            load(filesTemp{iProbe,iTrial});
            
            switch method
                case 'cbp'
                    dataProbe.data = [dataProbe.data,spikes.data];
                    dataProbe.dt   = 1 / spikes.Fs;
                case 'fmm'
                    spikesAll = psr_sst_spike_append(spikesAll,spikes);
                case 'kst'
                    dataProbe = [dataProbe,spikes.data]; %#ok
                    Fs = spikes.Fs;
                case 'ops'
                    spikesAll = psr_sst_spike_append(spikesAll,spikes);
                case 'ums'
                    spikesAll = psr_sst_spike_append(spikesAll,spikes);
            end
        end
        
        %% SPIKE SORT
        
        switch method
            case 'cbp' % WORK IN PROGRESS
                assignsAll = psr_sst_sorting_CBP(dataProbe,parameters);
                clear data_probe;
            case 'fmm'
                assignsAll = psr_sst_sorting_FMM(spikesAll,parameters);
            case 'kst'
                parameters.Fs = Fs;
                rez       = psr_sst_sorting_KST(dataProbe,parameters,savePath);
                spikesAll = psr_kst_convert2spikes(rez,dataProbe,parameters);
                spikesAll.Fs = rez.ops.fs;
                fclose('all');
                delete(rez.ops.fbinary);
                clear data_probe;
            case 'ops' % WORK IN PROGRESS
                assignsAll = psr_sst_sorting_OPS(spikesAll,parameters);
                clear spikesAll;
            case 'ums'
                assignsAll = psr_sst_sorting_UMS(spikesAll,parameters);
                clear spikesAll;
        end
        
        %% SPLIT DATA BACK INTO TRIALS
                
        nTrials = size(filesTemp,2);
        freq = []; % Initialize in case not loaded
        
        % Combine spikes structures if needed
        
        tf = strcmp(method,{'fmm','ums','ops'});
        tf = max(tf(:));
                
        if (tf)
            N = 0;
            spikesAll = [];
            for iTrial = 1:nTrials
                load(filesTemp{iProbe,iTrial});
                spikes.assigns = assignsAll(N(iTrial)+1:N(iTrial+1));
                spikesAll = psr_sst_spike_append(spikesAll,spikes);
                N = N + size(spikes.spiketimes,2);
            end
        end
        
        % Remove spikes based on RPVs
        spikesAll    = psr_sst_filter_rpv(spikesAll,parameters); 
        
        % Initialize arrays
        trials       = zeros(size(spikesAll.spiketimes),'int16');
        stimVoltsAll = zeros(nTrials,1);
        onsetTimes   = zeros(nTrials,1);
        offsetTime   = 0;
        offsetSpike  = 0;
        
        % Store data per trial
        for iTrial = 1:nTrials
            load(filesTemp{iProbe,iTrial}); % Load temporary trial data
            if (iTrial == 1); stimTimesAll = cell(nTrials,length(metadata.stimtimes)); end
            stimTimesAll(iTrial,:) = metadata.stimtimes;
            stimVoltsAll(iTrial)   = metadata.stimulus;
            if (~isempty(freq))
                freqAll(iTrial) = freq; %#ok
            end
            onsetTime   = offsetTime;
            offsetTime  = onsetTime + metadata.duration;
            onsetSpike  = offsetSpike + 1;
            offsetSpike = find(spikesAll.spiketimes <= offsetTime,1,'last');
            onsetTimes(iTrial) = onsetTime;
            trials(onsetSpike:offsetSpike) = iTrial;
        end
        
        % Merge clusters
        
        metadata.stimtimes  = stimTimesAll;
        metadata.trialonset = onsetTimes;
        metadata.stimulus   = stimVoltsAll;
        metadata.duration   = spikesAll.info.dur;
        
        spikesAll.trials = trials;
        spikesAll = psr_sst_cluster_merge(spikesAll,freqAll,metadata,parameters); % Merge clusters
        
        % Save
        
        filesSaved{iProbe} = saveFile(spikesAll,freqAll,metadata,parameters,savePath); % Save to MAT file
        save([savePath tempFileName],'filesSaved','-append'); % Update output file array
    end
    
    % Remove spikes based on high correlation across probes
    if (parameters.spikes.artifacts_removal)
        disp('Calculating correlations across probes...');
        psr_sst_artifact_correlation(filesSaved,filesTemp);
    end
    
    %% Delete temporary files
    delete([savePath tempFileName '.mat']);
    for iProbe = 1:nProbes
        for iTrial = 1:nTrials
            delete(filesTemp{iProbe,iTrial});
        end
    end
    
    disp(['Spike sorting completed. MAT file(s) saved to: "' savePath '"']);
    
end

end

function filename = saveFile(spikes,freq,metadata,parameters,savePath)

filename = [savePath ...
    'Spikes_'   metadata.subject                  ...
    '_'         strjoin(metadata.session, '-')    ...
    '_P'        num2str(metadata.probe,   '%02d') ...
    '_'         parameters.sorting.method];

if (length(metadata.stimulus) == 1)
    filename = [filename ...
        '_V'    num2str(metadata.stimulus, '%03d')];
end

freq       = orderfields(freq);       %#ok
parameters = orderfields(parameters); 
spikes     = orderfields(spikes);     %#ok

filename = [filename '.mat'];
save(filename,parameters.general.savelist{:});

end

function val = extractStringFromPath(filepath,token)
val      = 0;
filepath = fliplr(filepath);
k = strfind(filepath,'\');
k = k(2);
filepath = filepath(1:k);
filepath = fliplr(filepath);
elems = regexp(filepath, token, 'tokens', 'once');
if (~isempty(elems)); str = elems{1};
    if (~isempty(str)); val = str2double(elems{1}); end
end

end

%------------- END OF CODE --------------