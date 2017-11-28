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

if (~isfield(parameters.general,'configPath') || exist(parameters.general.configPath,'file') == 0)
    psr_parameter_default;
else                          
    run(parameters.general.configPath);
end
parameters = orderfields(parameters);

%% Constants
nElectrodes  = parameters.general.nelectrodes;
pattern      = parameters.general.filepattern;
ext          = parameters.general.extension;
precision    = 10^parameters.general.precision;
tempFileName = 'Temp_vars';

parameters.lfp.freqFields = ...
    {'artifacts', ...
    'cfg',        ...
    'dimord',     ...
    'freq',       ...
    'fsample',    ...
    'hdr',        ...
    'label',      ...
    'powspctrm',  ...
    'sampleinfo', ...
    'std',        ...
    'time',       ...
    'trialIDs',   ...
    'trial'};

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
    
    if (parameters.process.ap)
        
        stimTimes = cell(nTrials,2); % stimulus onset times
        k = psr_session_strcmp(parameters,'passive');
        loadPathActive  = loadPath(~k);
        loadPathPassive = loadPath( k);
        
        % Camera onset times
        
        % For active trials
        stimTimes(~k,1) = psr_cam_detection(loadPathActive);
        stimTimes(~k,2) = cellstr('interval');
        
        % Magnetic field artifacts [MFA]
        
        % For passive trials
        parameters.general.stims = stimuliConditions(k);
        stimTimes(k,1) = psr_mfa_detection(loadPathPassive,parameters);
        stimTimes(k,2) = cellstr('onset');
        
    end
    
    %% SPIKE + LFP DETECTION
    
    nProbes = zeros(nTrials,1);
    parameters.Fr = parameters.lfp.rsfactor * parameters.lfp.bp_upper;
    
    for iTrial = 1:nTrials
        
        % Check for existing files
        
        if (iTrial <= size(filesTemp,2))
            filesRaw = filesTemp(:,iTrial);
            k = ~cellfun(@isempty,filesRaw);
            filesRaw = filesRaw(k);
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
        ts_LFP           = [];
        
        data_mean = []; %%%% TEMP
        
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
                
                data_channel_raw = detrend(data_channel_raw, 'linear');
                data_channel_raw = psr_artifact_fft(data_channel_raw,parameters,Fs);
                
                %% SPIKING DATA
                
                if (parameters.process.spikes)
                    
                    %% Initialize arrays
                    if (iElectrode == 1)
                        ts_Spikes = [];
                        ts_Spikes.data = zeros(nElectrodes,length(data_channel_raw),'int16');
                    end
                    
                    % BAND-PASS FILTERING
                    cfg       = [];
                    cfg.Fs    = Fs;
                    cfg.order = parameters.spikes.bp_order;
                    cfg.lower = parameters.spikes.bp_lower;
                    cfg.upper = parameters.spikes.bp_upper;
                    data_channel = psr_bp_filter(data_channel_raw,cfg);
                    
                    ts_Spikes.data(iElectrode,:) = int16(precision * data_channel);
                    
                end
                
                %% LOCAL FIELD POTENTIAL
                if (parameters.process.lfp)
                    
                    % Resample raw signal
                    timestamps        = timestamps - timestamps(1);
                    [data,timestamps] = resample(data_channel_raw,timestamps,parameters.Fr); % Initial resample to avoid errors in ft_preprocessing
                    
                    if (iElectrode == 1)
                        ts_LFP.data = zeros(nElectrodes,length(data),'single');
                        ts_LFP.time = zeros(nElectrodes,length(data),'single');
                    end
                    
                    ts_LFP.data(iElectrode,:) = data;
                    ts_LFP.time(iElectrode,:) = timestamps;
                    
                end
                
                %%%% TEMP
                if (isempty(data_mean)); data_mean = data_channel_raw;
                else,                    data_mean = data_mean + data_channel_raw;
                end
                %%%%
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
            
            if (parameters.process.spikes); ts_Spikes.Fs = Fs; end % Hz, sampling rate of spike data
            
            % Save
            
            % Create metadata variable: to be expanded in future versions
            % with data from the electronics notebook
            
            metadata = [];
            metadata.subject   = subject;
            metadata.session   = session;
            metadata.dir       = loadPath;
            metadata.stimtimes = stimTimesTrial;
            metadata.duration  = sLength / Fs;
            metadata.stimulus  = stim;
            metadata.probe     = iProbe;
            
            % Save temporary MAT file
            filename = [savePath 'Temp_Probe_' num2str(iProbe,'%02d') '_Trial_' num2str(iTrial,'%02d') '.mat'];
            filesTemp{iProbe,iTrial} = filename;
            save(filename,'ts_Spikes','ts_LFP','metadata','parameters');
        end
        %%%% TEMP
        data_mean = data_mean / 64;
        save([savePath 'data_' num2str(iTrial,'%02d')],'data_mean','stimTimesTrial');
        %%%%
    end
    
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

else
    load([savePath tempFileName])
end

nProbes = size(filesTemp,1);
nTrials = size(filesTemp,2);

if (~exist('filesSaved','var')) % Initialize array if not loaded
    filesSaved = cell(nProbes,3); % Output files
end

%% Local field potential

if (parameters.process.lfp)
    
    for iProbe = 1:nProbes
        
        if (exist(filesSaved{iProbe,1},'file') && strcmp(filesSaved{iProbe,2},'LFP')); continue; end
        
        dataProbe = cell(nTrials,1);
        timeProbe = cell(nTrials,1);
        stimProbe = cell(nTrials,1);
        
        % Load all trials for probe
        for iTrial = 1:nTrials
            load(filesTemp{iProbe,iTrial},'ts_LFP','metadata');
            dataProbe{iTrial} = ts_LFP.data;
            timeProbe{iTrial} = ts_LFP.time;
            stimProbe{iTrial} = metadata.stimtimes;
        end
        
        load(filesTemp{iProbe,1},'parameters'); % Load parameters from first trial
        
        artifacts = [];
        artifacts = [artifacts,psr_lfp_artifact_detection_psd(dataProbe,parameters)]; %#ok
        artifacts = [artifacts,psr_lfp_artifact_detection_amp(dataProbe,parameters)]; %#ok
        [dataProbe,artifacts] = psr_lfp_artifact_removal(dataProbe,artifacts,parameters);
        
        % Initialize structure array
                
        freqFields = parameters.lfp.freqFields;
        nFields = length(freqFields);
        for iField = 1:nFields
            freqArray(nTrials).(freqFields{iField}) = []; %#ok
        end
        
        % Do LFP processing
        
        for iTrial = 1:nTrials
            freq = [];
            stimTimesTrial = stimProbe{iTrial};
            if (~isempty(stimTimesTrial))
                
                inputs            = [];
                inputs.data       = dataProbe{iTrial};
                inputs.timestamps = mean(timeProbe{iTrial},1);
                inputs.stimtimes  = stimTimesTrial;
                inputs.method     = 'tfa';
                
                load(filesTemp{iProbe,iTrial},'parameters');
                output = psr_lfp_wrapper(inputs,parameters);
                parameters     = output.parameters;
                freq           = output.freq;
                freq.artifacts = artifacts{iTrial};
            end
            
            % Save to structure array
            for iField = 1:nFields
                if isfield(freq,freqFields{iField})
                    freqArray(iTrial).(freqFields{iField}) = freq.(freqFields{iField});
                end
            end
        end
        
        % Save LFP output
        parameters.general.savelist = {'freq','parameters'};
        filesSaved{iProbe,1} = saveFile([],freqArray,metadata,parameters,savePath,true);
        filesSaved{iProbe,2} = 'LFP';
        save([savePath tempFileName],'filesSaved','-append'); % Update output file array
    end
end

%% Spike sorting

if (parameters.process.spikes)
    
    sortMethod = lower(parameters.sorting.method);
    
    for iProbe = 1:nProbes % Do clustering per probe across all stimulus conditions
        
        % Check if output file + temp file (after sorting) exist
        if (exist(filesSaved{iProbe},'file') && strcmp(filesSaved{iProbe,3},'SPK')); continue; end
        
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
        tic; % Measure duration of sorting
        
        % Initialize
        dataProbe = [];
        spikesAll = []; % structure to contain all 'spikes' structures from all trials
                
        %% MERGE DATA ACROSS TRIALS
        
        Fs_array = zeros(nTrials,1);
        
        for iTrial = 1:nTrials
            load(filesTemp{iProbe,iTrial});
            parameters.Fs = ts_Spikes.Fs;
            Fs_array(iTrial) = parameters.Fs;
            if (max(strcmp(sortMethod,{'cbp','kst','ops'}))) % Concatenate data
                dataProbe = [dataProbe,ts_Spikes.data]; %#ok
            end
            if (max(strcmp(sortMethod,{'fmm','iso','ops','ost','spc','ums'}))) % Spike detection
                spikes    = psr_sst_detection(ts_Spikes.data,parameters);
                spikesAll = psr_sst_spike_append(spikesAll,spikes);
            end
        end
        
        flag = sum(Fs_array ~= Fs_array(1)) ~= 0;
        if (flag); disp(['Sampling frequency mismatch. Skipping probe ' num2str(iProbe) '...']); continue;
        else; Fs = Fs_array(1);
        end
           
        %% SPIKE SORT
                
        switch sortMethod
            case 'cbp'
                psr_parameter_spikes; % TEMP
                input      = [];
                input.data = dataProbe; 
                input.dt   = 1 / Fs;
                clear dataProbe;
                assignsAll = psr_sst_sorting_CBP(input,parameters);
                spikesAll.assigns = assignsAll;
            case 'fmm'
                spikesAll.assigns = psr_sst_sorting_FMM(spikesAll,parameters);
            case 'iso'
                psr_parameter_spikes; % TEMP
                spikesAll.assigns = psr_sst_sorting_ISO(spikesAll,parameters);
            case 'kst'
                psr_parameter_default; % TEMP
                parameters.Fs = Fs;
                rez = psr_sst_sorting_KST(dataProbe,parameters,savePath);
                spikesAll = psr_kst_convert2spikes(rez,dataProbe,parameters);
                spikesAll.Fs = rez.ops.fs;
                fclose('all'); % Clean-up
                delete(rez.ops.fbinary);
                clear dataProbe;
            case 'ops' % WORK IN PROGRESS
                spikesAll.assigns = psr_sst_sorting_OPS(spikesAll,dataProbe,parameters);
            case 'ost'
                spikesAll.assigns = psr_sst_sorting_OST(spikesAll,parameters);
            case 'spc'
                spikesAll.assigns = psr_sst_sorting_SPC(spikesAll,parameters);
            case 'ums'
                spikesAll.assigns = psr_sst_sorting_UMS(spikesAll,parameters);
        end
        
        spikesAll.time = toc;
        
        %% SPLIT DATA BACK INTO TRIALS
        
        nTrials = size(filesTemp,2);
        freq = []; % Initialize in case not loaded
                
        % Initialize arrays
        trials       = zeros(size(spikesAll.spiketimes),'int16');
        stimVoltsAll = zeros(nTrials,1);
        onsetTimes   = zeros(nTrials,1);
        offsetTime   = 0;
        offsetSpike  = 0;
        
        % Store data per trial
        for iTrial = 1:nTrials
            load(filesTemp{iProbe,iTrial},'metadata'); % Load temporary trial data
            if (iTrial == 1); stimTimesAll = cell(nTrials,length(metadata.stimtimes)); end
            stimTimesAll(iTrial,:) = metadata.stimtimes;
            stimVoltsAll(iTrial)   = metadata.stimulus;
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
        spikesAll = psr_sst_cluster_merge(spikesAll,parameters); % Merge clusters
        
        % Save
        
        if (exist(filesSaved{iProbe},'file') > 0) % Load LFP parameters
            parametersTemp = parameters;
            load(filesSaved{iProbe},'parameters');
            parametersTemp.lfp = parameters.lfp;
            parameters = parametersTemp;
        end
        
        metadata.sessionIndex = parameters.general.sessionIndex;
        parameters.general.savelist = {'spikes','metadata','parameters'}; % What variables to save in output MAT file
        filesSaved{iProbe,1} = saveFile(spikesAll,[],metadata,parameters,savePath,true); % Save to MAT file
        filesSaved{iProbe,3} = 'SPK';
        save([savePath tempFileName],'filesSaved','-append'); % Update output file array
    end
       
    %% Stability check (for development)
    if (parameters.develop.comparison)
        [filesSaved,filesTemp] = psr_stability_check(filesSaved,filesTemp); %#ok
        save([savePath tempFileName],'filesSaved','filesTemp','-append'); % Update output file array
        return; % Do another comparison
    end
    
    %% Remove spikes based on high correlation across probes
    if (parameters.spikes.artifacts_removal)
        disp('Calculating correlations across probes...');
        psr_sst_artifact_correlation(filesSaved,filesTemp);
    end
    
    %% Calculate cluster features
    for iProbe = 1:nProbes
        disp(['Calculating cluster features for probe ' num2str(iProbe) '...']);
        filename = filesSaved{iProbe};
        load(filename);
        if (~isfield(spikes,'clusters'))
            spikes.clusters = psr_sst_cluster_features(spikes,freq,metadata,parameters);
        end
        save(filename,'spikes','-append');
    end
    
    disp(['Spike sorting completed. MAT file(s) saved to: "' savePath '"']);
    
end

%% Delete temporary files
delete([savePath tempFileName '.mat']);
for iProbe = 1:nProbes
    for iTrial = 1:nTrials
        delete(filesTemp{iProbe,iTrial});
    end
end

end

function filename = saveFile(spikes,freq,metadata,parameters,savePath,append)

if (nargin < 6); append = false; end

saveList = parameters.general.savelist;

filename = [savePath ...
    'Spikes_'   metadata.subject                  ...
    '_'         strjoin(metadata.session, '-')    ...
    '_P'        num2str(metadata.probe,   '%02d') ...
    '_'         upper(parameters.sorting.method)];

if (~isempty(freq));   freq   = orderfields(freq);   end %#ok
if (~isempty(spikes)); spikes = orderfields(spikes); end %#ok
parameters = orderfields(parameters); %#ok

filename = [filename '.mat'];
if (append && exist(filename,'file') > 0); save(filename,saveList{:},'-append');
else,                                      save(filename,saveList{:});
end

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