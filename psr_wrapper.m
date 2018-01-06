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
        loadPath{iTrial} = [loadPath{iTrial}, '\']; %#ok
    end
end

if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

%% Convert initial parameters
temp.general = parameters;
parameters   = temp; clear temp;

%% Load processing parameters
parameters = loadParameters(parameters);

%% Constants
nChans       = parameters.general.nelectrodes;
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
    iProbe = extractStringFromPath(filename,'Probe_(\d*)_');
    iTrial = extractStringFromPath(filename,'Trial_(\d*).');
    filesTemp{iProbe,iTrial} = filename;
end
filesTempFlag = cellfun(@isempty,filesTemp);

fileTemp = dir([savePath tempFileName '.mat']);
fileTemp = char(fileTemp.name);
fileTempFlag = isempty(fileTemp);

%% Start data processing

if (isempty(filesTempFlag) || max(filesTempFlag(:)))
    
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
    
    stimTimes = cell(nTrials,2); % stimulus onset times
    
    % Different detection methods for stimulus onset times. This is heavily
    % dependent on the experimental paradigm.
    
    %% ACTIVE VS. PASSIVE
    
    if (parameters.experimental.ap.process)
        
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
    parameters.Fr = parameters.lfp.rsfactor * parameters.lfp.bp.upper;
    
    for iTrial = 1:nTrials
        
        % Initialize
        
        filesRaw        = filesRawAll{iTrial};
        nFiles          = length(filesRaw);
        nProbes(iTrial) = nFiles / parameters.general.nelectrodes;
        stimTimesTrial  = stimTimes(iTrial,:);
        stim            = stimuliConditions(iTrial);
        ts_Spikes       = [];
        ts_LFP          = [];
        
        for iProbe = 1:nProbes(iTrial)
            
            % Check if TEMP file already exists
            if (~isempty(filesTempFlag) && ...
                    iProbe <= size(filesTempFlag,1) && ...
                    iTrial <= size(filesTempFlag,2))
                if (~filesTempFlag(iProbe,iTrial)); continue; end
            end
            
            for iChan = 1:nChans
                
                % Load Open-Ephys data
                
                iFile = (iProbe - 1) * nChans + iChan;
                
                filename = [loadPath{iTrial} filesRaw{iFile}];
                filename = strtrim(filename);
                
                try % Load CONTINUOUS files [microvolts]
                    disp(['Loading ' filename '...']);
                    [dataChanRaw, timestamps, info] = load_open_ephys_data_faster(filename);
                    timestamps = timestamps / info.header.sampleRate;
                catch
                    [dataChanRaw, timestamps, info] = load_open_ephys_data(filename);
                end
                
                % Store variables
                if (iChan == 1) % Initialize
                    Fs_array = zeros(nChans,1);
                    Ls_array = zeros(nChans,1);
                end
                
                Fs = info.header.sampleRate;
                Fs_array(iChan) = Fs; % Sampling frequencies should be equal
                Ls_array(iChan) = length(dataChanRaw);
                
                dataChanRaw = detrend(dataChanRaw, 'linear');
                
                if (iChan == 1); dataChanRawSum = dataChanRaw;
                else,            dataChanRawSum = dataChanRaw + dataChanRawSum;
                end
                
                dataChanRaw = psr_artifact_fft(dataChanRaw,parameters,Fs);
                
                if (parameters.experimental.ap.process ...
                        && parameters.experimental.ap.diff && stim > 0)
                    dataChanRaw = psr_artifact_diff(dataChanRaw,parameters,Fs);
                end
                
                %% SPIKING DATA
                
                if (parameters.process.spikes)
                    
                    % BAND-PASS FILTERING
                    cfg       = [];
                    cfg.Fs    = Fs;
                    cfg.order = parameters.spikes.bp.order;
                    cfg.lower = parameters.spikes.bp.lower;
                    cfg.upper = parameters.spikes.bp.upper;
                    
                    dataChanFiltered = psr_bp_filter(dataChanRaw,cfg);
                    
                    if (iChan == 1)
                        ts_Spikes.data = zeros(nChans,length(dataChanFiltered),'int16');
                    end
                    
                    ts_Spikes.data(iChan,:) = int16(precision * dataChanFiltered);
                    
                end
                
                %% LOCAL FIELD POTENTIAL
                
                if (parameters.process.lfp)
                    
                    % Resample raw signal
                    timestamps        = timestamps - timestamps(1);
                    [data,timestamps] = resample(dataChanRaw,timestamps,parameters.Fr); % Initial resample to avoid errors in ft_preprocessing
                    
                    if (iChan == 1)
                        ts_LFP.data = zeros(nChans,length(data),'single');
                        ts_LFP.time = zeros(nChans,length(data),'single');
                    end
                    
                    ts_LFP.data(iChan,:) = data;
                    ts_LFP.time(iChan,:) = timestamps;
                    
                end
                
            end
            
            clear data dataChanRaw
            
            % Check data
            
            flag = sum(Fs_array ~= Fs_array(1)) ~= 0;
            if (flag); disp('Sampling frequency mismatch'); continue;
            else; Fs = Fs_array(1);
            end
            
            flag = sum(Ls_array ~= Ls_array(1)) ~= 0;
            if (flag); disp('Data length mismatch'); continue;
            else; sLength = Ls_array(1);
            end
            
            if (parameters.process.spikes)
                ts_Spikes.Fs = Fs; % Hz, sampling rate of spike data
            end
            
            % Save
            
            % Create metadata variable: to be expanded in future versions
            % with data from the electronics notebook
            
            metadata = [];
            metadata.subject    = subject;
            metadata.session    = session;
            metadata.stimtimes  = stimTimesTrial;
            metadata.duration   = sLength / Fs;
            metadata.stimulus   = stim;
            metadata.probe      = iProbe;
            
            if (parameters.experimental.ap.process ...
                    && parameters.experimental.ap.offset && stim > 0)
                metadata.stimoffset = psr_sst_artifact_stimoffset(dataChanRawSum,stimTimesTrial,Fs);
            end
            
            % Save temporary MAT file
            filename = [savePath 'Temp_Probe_' num2str(iProbe,'%02d') '_Trial_' num2str(iTrial,'%02d') '.mat'];
            filesTemp{iProbe,iTrial} = filename;
            save(filename,'ts_Spikes','ts_LFP','metadata','parameters');
        end
    end
    
end

% Temporarily save all workspace variables. Useful when continuing after early stopping

if (fileTempFlag)
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
    save([savePath tempFileName]);
else
    load([savePath tempFileName])
end

nProbes = size(filesTemp,1);
nTrials = size(filesTemp,2);

% Sort trials

sessions  = parameters.general.sessionIndex;
nSessions = length(unique(sessions));
trialIDs  = 1:nTrials;
itr       = 1;
for iSession = 1:nSessions
    stimSession = stimuliConditions(sessions == iSession);
    [~,I] = sort(stimSession);
    n = length(stimSession);
    id = trialIDs(itr:itr+n-1);
    trialIDs(itr:itr+n-1) = id(I);
    itr = itr + n;
end

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
                parameters = loadParameters(parameters); % TEMP
                output = psr_lfp_wrapper(inputs,parameters);
                parameters     = output.parameters;
                freq           = output.freq;
                freq.artifacts = artifacts{iTrial};
            end
            
            % Save to structure array
            for iField = 1:nFields
                if isfield(freq,freqFields{iField})
                    freqArray(trialIDs(iTrial)).(freqFields{iField}) = freq.(freqFields{iField});
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

parameters = loadParameters(parameters); % Reload parameters

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
        
        disp(['Spike sorting probe ' num2str(iProbe,'%02d') '...']);
        tSort = tic; % Measure duration of spike detection and sorting
        
        % Initialize
        dataProbe = [];
        spikesAll = []; % structure to contain all 'spikes' structures from all trials
        
        %% MERGE DATA ACROSS TRIALS
        
        Fs_array = zeros(nTrials,1);
        
        for iTrial = 1:nTrials
            trialID = trialIDs(iTrial);
            load(filesTemp{iProbe,trialID});
            parameters.Fs = ts_Spikes.Fs;
            Fs_array(iTrial) = parameters.Fs;
            
            if (parameters.develop.timing) % For development
                parameters = loadParameters(parameters);
                n = round(parameters.Fs * parameters.develop.time);
                ts_Spikes.data = ts_Spikes.data(:,1:n);
            end
            
            if (max(strcmp(sortMethod,{'cbp','kst','ops'}))) % Concatenate data
                dataProbe = [dataProbe,ts_Spikes.data]; %#ok
                SPIKE_DETECTION = false;
            end
            
            if (max(strcmp(sortMethod,{'fmm','iso','kfm','ops','ost','spc','ums'}))) % Spike detection
                spikes    = psr_sst_detection(ts_Spikes.data,parameters);
                spikesAll = psr_sst_spike_append(spikesAll,spikes);
                SPIKE_DETECTION = true;
            end
        end
        
        flag = sum(Fs_array ~= Fs_array(1)) ~= 0;
        if (flag); disp(['Sampling frequency mismatch. Skipping probe ' num2str(iProbe) '...']); continue;
        else; parameters.Fs = Fs_array(1);
        end
        
        %% SPIKE SORT
        
        parameters = loadParameters(parameters); % re-load, in case changes were made
        
        switch sortMethod
            case 'cbp'; spikesAll         = psr_sst_sorting_CBP(dataProbe,parameters);
            case 'fmm'; spikesAll.assigns = psr_sst_sorting_FMM(spikesAll,parameters);
            case 'iso'; spikesAll.assigns = psr_sst_sorting_ISO(spikesAll,parameters);
            case 'kfm'; spikesAll.assigns = psr_sst_sorting_KFM(spikesAll,parameters);
            case 'kst'; spikesAll         = psr_sst_sorting_KST(dataProbe,parameters,savePath);
            case 'ops'; spikesAll.assigns = psr_sst_sorting_OPS(spikesAll,dataProbe,parameters);
            case 'ost'; spikesAll.assigns = psr_sst_sorting_OST(spikesAll,parameters);
            case 'spc'; spikesAll.assigns = psr_sst_sorting_SPC(spikesAll,parameters);
            case 'ums'; spikesAll.assigns = psr_sst_sorting_UMS(spikesAll,parameters);
        end
        
        if ~isempty(spikesAll) || parameters.develop.timing
            
            if (parameters.develop.timing) % For development
                tEnd = toc(tSort); % Spike sorting time
                filename = [savePath 'timings.mat'];
                if (exist(filename,'file'))
                    load(filename);
                else
                    timings = [];
                    methods = [];
                end
                methods{end+1} = sortMethod;
                timings(end+1) = tEnd;
                save([savePath 'timings.mat'],'methods','timings');
                continue;
            end
            
            %% SET TRIAL INDEX FOR EACH SPIKE
            
            spikesAll.Fs = parameters.Fs; % Store sampling frequency
            
            % Initialize arrays
            
            nTrials     = size(filesTemp,2);
            trials      = zeros(size(spikesAll.spiketimes),'int16');
            offsetTime  = 0;
            offsetSpike = 0;
            
            for iTrial = 1:nTrials
                trialID = trialIDs(iTrial);
                load(filesTemp{iProbe,trialID},'metadata','parameters'); % Load temporary trial data
                onsetTime   = offsetTime;
                offsetTime  = onsetTime + metadata.duration;
                onsetSpike  = offsetSpike + 1;
                offsetSpike = find(spikesAll.spiketimes <= offsetTime,1,'last');
                trials(onsetSpike:offsetSpike) = iTrial;
            end
            
            spikesAll.trials = trials;
                        
            %% SAVE
            
            if (exist(filesSaved{iProbe},'file') > 0) % Load LFP parameters
                parametersTemp = parameters;
                load(filesSaved{iProbe},'parameters');
                parametersTemp.lfp = parameters.lfp;
                parameters = parametersTemp;
            end
            
            spikesAll.info.detected = SPIKE_DETECTION;
            
            load(filesTemp{iProbe,1},'metadata'); % Load metadata from first trial
            parameters.general.savelist = {'spikes','parameters'}; % What variables to save in output MAT file
            filesSaved{iProbe,1} = saveFile(spikesAll,[],metadata,parameters,savePath,true); % Save to MAT file
            
        end
        
        filesSaved{iProbe,3} = 'SPK';
        save([savePath tempFileName],'filesSaved','-append'); % Update output file array
    end
    
    if (parameters.develop.timing); return; end
    
    %% Calculate some metrics if needed
    psr_sst_spikes_info(filesSaved,filesTemp);
    
    %% Calculate correlation of spikes across probes
    disp('Calculating correlations across probes...');
    psr_sst_corr_global(filesSaved,filesTemp);
    
    %% Calculate cluster features
    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    for iProbe = 1:nProbes
        disp(['Calculating cluster features for probe ' num2str(iProbe,'%02d') '...']);
        spikes = [];
        freq   = [];
        filename = filesSaved{iProbe};
        load(filename);
        if (~isempty(spikes))
            if (~isfield(spikes,          'delete')); spikes = psr_sst_filter_spikes(spikes,parameters); end % Filter spikes
            if (~isfield(spikes,   'assigns_prior')); spikes = psr_sst_cluster_merge(spikes,parameters); end % Merge clusters
            if (~isfield(spikes.clusters,'metrics')); spikes.clusters.metrics = psr_sst_cluster_features(spikes,freq,parameters); end % Calculate cluster metrics
        end
        save(filename,'spikes','-append');
    end
    warning('on', MSGID);
    
    disp(['Spike sorting completed. MAT file(s) saved to: "' savePath '"']);
    
end

%% Merge metadata across trials

% Initialize arrays
nProbes = size(filesTemp,1);
nTrials = size(filesTemp,2);

for iProbe = 1:nProbes
    
    stimVoltsAll = zeros(nTrials,1);
    onsetTimes   = zeros(nTrials,1);
    offsetTime   = 0;
    stimOffsets  = [];
    
    for iTrial = 1:nTrials
        
        trialID = trialIDs(iTrial);
        load(filesTemp{iProbe,trialID},'metadata','parameters'); % Load temporary trial data
        
        if (iTrial == 1); stimTimesAll = cell(nTrials,length(metadata.stimtimes)); end
        stimTimesAll(iTrial,:) = metadata.stimtimes;
        stimVoltsAll(iTrial)   = metadata.stimulus;
        if (isfield(metadata,'stimoffset')); stimOffsets(iTrial,:) = sort(metadata.stimoffset)'; end
        
        onsetTime  = offsetTime;
        offsetTime = onsetTime + metadata.duration;
        onsetTimes(iTrial) = onsetTime;
        
    end
    
    metadata.stimtimes  = stimTimesAll;
    metadata.trialonset = onsetTimes;
    metadata.stimulus   = stimVoltsAll;
    metadata.duration   = offsetTime;
    if (~isempty(stimOffsets)); metadata.stimoffset = stimOffsets; end
    
    parameters.general.savelist = {'metadata'}; % What variables to save in output MAT file
    saveFile([],[],metadata,parameters,savePath,true); % Save to MAT file
    
end

%% Stability check (for development)
if (parameters.develop.comparison)
    filesSaved = psr_stability_check(filesSaved,filesTemp,sortMethod); %#ok
    save([savePath tempFileName],'filesSaved','-append'); % Update output file array
    return;
end

%% Delete temporary files
delete([savePath tempFileName '.mat']);
for iProbe = 1:nProbes
    for iTrial = 1:nTrials
        delete(filesTemp{iProbe,iTrial});
    end
end

end

function parameters = loadParameters(parameters)

if (~isfield(parameters.general,'configPath') || exist(parameters.general.configPath,'file') == 0)
    psr_parameter_default;
else
    run(parameters.general.configPath);
end

parameters = orderfields(parameters);

end

function filename = saveFile(spikes,freq,metadata,parameters,savePath,append)

if (nargin < 6); append = false; end

saveList = parameters.general.savelist;

filename = [savePath ...
    'Spikes_'   metadata.subject                  ...
    '_'         strjoin(metadata.session, '-')    ...
    '_P'        num2str(metadata.probe,   '%02d')];

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