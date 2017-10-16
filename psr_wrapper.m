function psr_wrapper(parameters) 

% PSR_WRAPPER - Wrapper function that loads raw extracellular data,
% performs all data processing steps and saves the results to a MAT file.
% The processing pipeline consists of spike, local field potential and
% stimulus detection, as well as spike sorting.
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

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% email address: t.s.n.brouns@gmail.com
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

if (~isfield(parameters,'loadPathSub')); loadPath = uigetdir([],'Select data folder'); % Open folder selection dialogue
else,                                    loadPath = parameters.loadPathSub;
end

if (~iscell(loadPath)); loadPath = {loadPath}; end

numtrials = length(loadPath);
for iTrial = 1:numtrials
    if (loadPath{iTrial}(end) ~= '\' && ~isempty(loadPath{iTrial}))
        loadPath{iTrial} = [loadPath{iTrial}, '\'];
    end
end

if (~isfield(parameters,'savePathSub')); savePath = [];
else,                                    savePath = parameters.savePathSub;
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

nelectrodes  = parameters.general.nelectrodes;
pattern      = parameters.general.filepattern;
ext          = parameters.general.extension;
method       = lower(parameters.sorting.method);
tempFileName = 'Temp_vars';

%% Check if chosen spike sorting method requires spike detection

tf = strcmp(method,{'fmm','ums','ops'});
tf = max(tf(:));
if (tf); SPIKE_DETECTION = 1;
else;    SPIKE_DETECTION = 0;
end

tf = strcmp(method,{'kst','cbp','ops'});
tf = max(tf(:));
if (tf); SAVE_RAW_DATA = 1;
else;    SAVE_RAW_DATA = 0;
end

%% Check if TEMP files exist in save folder

files_all  = cell(0,0);
files_temp = dir([savePath 'Temp_Probe_*_Trial_*.mat']);
files_temp = char(files_temp.name);
nfiles     = size(files_temp,1);
for iFile = 1:nfiles
    filename = [savePath files_temp(iFile,:)];
    probe    = extractStringFromPath(filename,'Probe_(\d*)_');
    trial    = extractStringFromPath(filename,'Trial_(\d*).');
    files_all{probe,trial} = filename;
end
file_temp = dir([savePath tempFileName '.mat']);
file_temp = char(file_temp.name);

%% Start data processing

if (isempty(files_all) || isempty(file_temp))
    
    %% Create metadata variable: to be expanded in future versions with data
    % from the electronics notebook
    
    metadata = [];
    metadata.subject  = subject;
    metadata.session  = session;
    metadata.dir      = loadPath;
    clear subject session
    
    %% Extract stimulus conditions from folder name
    
    stimuliConditions = zeros(numtrials,1);
    
    for iTrial = 1:numtrials
        filepath = loadPath{iTrial};
        stimuliConditions(iTrial) = extractStringFromPath(filepath,['_(\d*)' parameters.general.trialpattern]);
    end
    
    %% Find files
    
    files_unsorted = cell(numtrials,1);
    for iTrial = 1:numtrials
        files_struct = dir([loadPath{iTrial} '\*' pattern '*' ext]);
        if (size(files_struct,1) == 0); return; end
        files_unsorted{iTrial} = char(files_struct.name);
    end
    
    %% Sort raw data files in correct chronological order 
    
    files_trials = cell(numtrials,1);
    
    for iTrial = 1:numtrials
        numfiles = length(files_unsorted{iTrial}(:,1));
        files_trial = cell(numfiles,1);
        for iFile = 1:numfiles
            filename   = files_unsorted{iTrial}(iFile,:);
            filename   = strtrim(filename);
            k          = strfind(filename,pattern) + length(pattern);
            [~,name,~] = fileparts(filename);
            id         = str2double(name(k:end));
            files_trial{id} = filename;
        end
        files_trial = files_trial(~cellfun('isempty',files_trial)); % remove empty cells
        files_trials{iTrial} = files_trial;
    end
    
    %% STIMULUS ONSET DETECTION (MAGNETIC FIELD ARTIFACTS [MFA])
    
    MFAtimes = cell(numtrials,2); % stimulus onset times
    
    if (parameters.process.mfa)
        
        % Files to load
        pattern     = 'ADC';
        fileIndices = [4 7];
        
        for iTrial = 1:numtrials
            
            % Load new files
            
            files = dir([loadPath{iTrial} '\*' pattern '*' ext]);
            files = char(files.name);
            
            if (size(files,1) >= max(fileIndices))
                files  = files(fileIndices,:);
                files  = strtrim(files);
                nfiles = length(fileIndices);
            elseif (size(files,1) > 0)
                files  = files(1,:);
                nfiles = 1;
            else
                continue;
            end
            
            % Initialize arrays
            
            data_channels =  cell(nfiles,1);
            Fs_array      = zeros(nfiles,1);
            
            % Filter raw data
            
            for iFile = 1:nfiles 
                
                file = [loadPath{iTrial} files(iFile,:)]; % Filename
                [data_ADC, ~, info] = load_open_ephys_data(file); % Load CONTINUOUS files [microvolts]
                
                Fs = info.header.sampleRate; % Sampling rate in Hz
                Fs_array(iFile) = Fs; % Save in array, because we have multiple sampling frequencies (should all be equal)
                
                % Band-pass filter
                
                [B,A] = butter(parameters.mfa.bp_order,[parameters.mfa.bp_lower parameters.mfa.bp_upper]/(Fs/2),'bandpass');
                data_ADC = filtfilt(B,A,data_ADC); % Zero-phase digital filtering
                
                % Normalize
                
                data_ADC = data_ADC / max(data_ADC);
                data_channels{iFile} = single(data_ADC); % convert to single type
                
            end
            
            parameters.Fs = mode(Fs_array); % Just take most frequently occuring sampling frequency
            
            if (stimuliConditions(iTrial) > 0) % If we have a stimulus
                [MFAtimes{iTrial,1},MFAtimes{iTrial,2}] = psr_mfa_detection(data_channels,parameters);
            else % Randomly select MFA times in order to get non-stimulus control data
                dur = length(data_channels{1}) / Fs;
                MFAtimes{iTrial,1} = sort(dur * rand(1,parameters.mfa.control));
                MFAtimes{iTrial,2} = sort(dur * rand(1,parameters.mfa.control));
            end
        end
        
    end
    
    %% SPIKE + LFP DETECTION
    
    %% Per probe
    
    nprobes = zeros(numtrials,1);
    parameters.Fr = parameters.lfp.rsfactor * parameters.lfp.bp_upper;
    
    for iTrial = 1:numtrials
        
        % Check for existing files
        
        if (iTrial <= size(files_all,2))
            files_trial = cell2mat(files_all(:,iTrial));
            numfiles    = size(files_trial,1);
            nprobes     = size(files_trials{iTrial},1) / parameters.general.nelectrodes;
            if (numfiles == nprobes); continue; end
        end
        
        % Initialize
        
        files_trial     = files_trials{iTrial};
        numfiles        = length(files_trial);
        nprobes(iTrial) = numfiles / parameters.general.nelectrodes;
        MFAtimesTrial   = MFAtimes(iTrial,:);
        stim            = stimuliConditions(iTrial);
        
        for iProbe = 1:nprobes(iTrial)
            
            for iElectrode = 1:nelectrodes
                
                % Load Open-Ephys data
                
                iFile = (iProbe - 1) * nelectrodes + iElectrode;
                
                filename = [loadPath{iTrial} files_trial{iFile}];
                filename = strtrim(filename);
                
                [data_channel_raw, timestamps, info] = load_open_ephys_data(filename); % data in microvolts
                
                % Store sampling frequency
                Fs = info.header.sampleRate;
                if (iElectrode == 1); Fs_array = zeros(nelectrodes,1); end
                Fs_array(iElectrode) = Fs; % Sampling frequencies should be equal
                
                data_channel_raw = psr_artifact_fft(data_channel_raw,parameters,Fs);
                
                %%%% TEMP
%                 % FFT
%                 Y = fft(data_channel_raw);
%                 L = length(data_channel_raw);
%                 P2 = abs(Y/L);
%                 P1 = P2(1:L/2+1);
%                 P1(2:end-1) = 2*P1(2:end-1);
%                 f = Fs*(0:(L/2))/L;
%                 figure(fig); plot(f,P1);
%                 xlabel('Frequency [Hz]')
%                 ylabel('Amplitude');
%                 ylim([0 3]);
%                 export_fig([savePath 'FFT_' num2str(iTrial) '_' num2str(iProbe) '_' num2str(iElectrode)]);
                %%%%
                
                if (parameters.process.spikes)
                    
                    %% Initialize arrays
                    if (iElectrode == 1)
                        data_probe_sst = zeros(nelectrodes,length(data_channel_raw),'single');
                        if (iProbe == 1 && parameters.spikes.artifacts_removal)
                            data_trial_sst = zeros(nprobes(iTrial)*nelectrodes,length(data_channel_raw),'single');
                        end
                    end
                    
                    %% BAND-PASS FILTERING SPIKE SORTING
                    
                    [B,A]        = butter(parameters.spikes.bp_order,[parameters.spikes.bp_lower parameters.spikes.bp_upper]/(Fs/2),'bandpass');
                    data_channel = filtfilt(B,A,data_channel_raw); % Zero-phase digital filtering
                    data_channel = single(data_channel); % convert to single type
                    
                    data_probe_sst(iElectrode,:) = detrend(data_channel, 'linear');
                    
                end
                
                %% LOCAL FIELD POTENTIAL
                
                % Resample raw signal
                                
                timestamps        = timestamps - timestamps(1);
                [data,timestamps] = resample(data_channel_raw,timestamps,parameters.Fr); % Initial resample to avoid errors in ft_preprocessing
                
                if (iElectrode == 1)
                    data_probe_lfp = zeros(nelectrodes,length(data),'single');
                    time_probe_lfp = zeros(nelectrodes,length(data),'single');
                end
                
                time_probe_lfp(iElectrode,:) = timestamps;
                data_probe_lfp(iElectrode,:) = data;
            end
            
            %%%% TEMP
            %             figure(fig);
            %             plot((1:size(data_probe_sst,2))/Fs,data_probe_sst);
            %             xlabel('Time [s]');
            %             ylabel('Voltage');
            %             savefig([savePath 'Signal_' num2str(iTrial) '_' num2str(iProbe) '.fig']);
            %%%%
            
            clear data_channel data_channel_raw
            
            flag = sum(Fs_array ~= Fs_array(1)) ~= 0;
            if (flag); disp('Sampling frequency mismatch'); continue;
            else; Fs = Fs_array(1);
            end
            
            spikes = [];
            spikes.params.Fs = Fs;  % Hz, sampling rate of spike data
            
            if (parameters.process.spikes)
                
                % SPIKE DETECTION
                
                nsamples         = 60 * parameters.spikes.twin * Fs; % cut data in sections
                nsamples_total   = size(data_probe_sst,2);
                nsection         =  ceil(nsamples_total / nsamples);
                nsamples_section = floor(nsamples_total / nsection); % process data in sections
                
                if (SPIKE_DETECTION)
                    iStart = 1;
                    for iSection = 1:nsection
                        data_section = data_probe_sst(:,iStart:iStart + nsamples_section - 1);
                        data_section = {data_section'};
                        if (stim > 0 && parameters.spikes.artifacts_combine) % magnetic field artifact detection
                            spikes = ss_mfa_detection_raw(data_section{1},spikes,iStart);
                        end
                        spikes = psr_sst_detect(data_section,spikes);
                        iStart = iStart + nsamples_section;
                    end
                end
                
                if (SAVE_RAW_DATA)
                    spikes.data = data_probe_sst;
                end
                
                if (parameters.spikes.artifacts_removal)
                    data_trial_sst(nelectrodes*(iProbe-1)+1:nelectrodes*iProbe,:) = data_probe_sst;
                end
            end
            
            spikes.stimtimes = MFAtimesTrial;
            
            % LOCAL FIELD POTENTIAL
            
            freq = [];
            if (parameters.process.lfp && ~isempty(spikes.stimtimes{1})) % TODO: REMOVE STIMULUS ONSET TIME CONDITION 
                timestamps = mean(time_probe_lfp,1);
                stimtimes  = MFAtimes{iTrial,1};
                [data,artifacts] = psr_artifact_removal_lfp(data_probe_lfp,parameters);
                [freq,parameters] = psr_lfp_wrapper(data,timestamps,stimtimes,parameters);
                freq.artifacts = artifacts;
            end
            
            % Save
            
            metadata.stimulus = stim;
            metadata.probe    = iProbe;
            
            if (parameters.process.spikes)
                % Save temporary MAT file
                filename = [savePath 'Temp_Probe_' num2str(iProbe) '_Trial_' num2str(iTrial) '.mat'];
                files_all{iProbe,iTrial} = filename;
                save(filename,'spikes','metadata','freq','parameters');
            else
                % Save LFP output
                saveFile(spikes,[],metadata,freq,parameters,savePath);
            end
        end
        
        %% Filter spiking data across all probes
                        
        if (parameters.process.spikes)
            
            filesTrial = files_all(:,iTrial);
            
            if (stim > 0 && parameters.spikes.artifacts_combine)
                MFAtimesTrial = ss_mfa_combine(filesTrial,MFAtimesTrial); % Magnetic field artifact combination
                psr_mfa_save(filesTrial,MFAtimesTrial)
            end
            
            %% Global median subtraction
            if (parameters.spikes.artifacts_subtract)
                psr_sst_artifact_subtraction(filesTrial, data_trial_sst);
            end
            
            %% Artifact removal
            if (SPIKE_DETECTION && parameters.spikes.artifacts_removal)
                psr_sst_artifact_removal(filesTrial, data_trial_sst, MFAtimesTrial);
            end
        end
    end
    
    if (parameters.process.spikes)
        keepvars = {...
            'parameters',        ...
            'files_all',         ...
            'files_trials',      ...
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

if (parameters.process.spikes)
    
    method  = lower(parameters.sorting.method);
    nprobes   = size(files_all,1);
    ntrials = size(files_all,2);
    
    for iProbe = 1:nprobes % Do clustering per probe across all stimulus conditions
        
        % Check data
        
        MISSING_FILE = false;
        for iTrial = 1:ntrials
            if (~exist(files_all{iProbe,iTrial},'file')); MISSING_FILE = true; break; end
        end
        if (MISSING_FILE) % If not all files are present
            for iTrial = 1:ntrials
                filename = files_all{iProbe,iTrial};
                if (exist(filename,'file')); delete(filename); end % Delete temporary spikes file
            end
            disp(['Skipping probe ' num2str(iProbe) '...']);
            continue; % Move on to next probe
        end
        
        disp(['Spike sorting probe ' num2str(iProbe) '...']);
        
        % Initialize
        
        switch method
            case 'cbp'
                data_probe      = []; %
                data_probe.data = [];
            case {'fmm','ops','ums'}
                spikesAll = []; % structure to contain all 'spikes' structures from all trials
            case 'kst'
                data_probe = [];
        end
        
        %% MERGE DATA ACROSS TRIALS
        
        for iTrial = 1:ntrials
            
            load(files_all{iProbe,iTrial});
            
            switch method
                case 'cbp'
                    data_probe.data = [data_probe.data,spikes.data];
                    data_probe.dt   = 1 / spikes.params.Fs;
                case 'fmm'
                    spikesAll = psr_sst_spike_append(spikesAll,spikes);
                case 'kst'
                    data_probe = [data_probe,spikes.data]; %#ok
                    Fs = spikes.params.Fs;
                case 'ops'
                    spikesAll = psr_sst_spike_append(spikesAll,spikes);
                case 'ums'
                    spikesAll = psr_sst_spike_append(spikesAll,spikes);
            end
        end
        
        data_probe = real(data_probe); % TEMP
        
        %% SPIKE SORT
        
        switch method
            case 'cbp' % WORK IN PROGRESS
                assigns_all = psr_sst_sorting_CBP(data_probe,parameters);
                clear data_probe;
            case 'fmm'
                assigns_all = psr_sst_sorting_FMM(spikesAll,parameters);
            case 'kst'
                parameters.Fs = Fs;
                rez = psr_sst_sorting_KST(data_probe,parameters,savePath);
                spikesAll = psr_kst_convert2spikes(rez,data_probe,parameters);
                spikesAll.params.Fs = rez.ops.fs;
                fclose('all');
                delete(rez.ops.fbinary);
                clear data_probe;
            case 'ops' % WORK IN PROGRESS
                assigns_all = psr_sst_sorting_OPS(spikesAll,parameters);
                clear spikesAll;
            case 'ums'
                assigns_all = psr_sst_sorting_UMS(spikesAll,parameters);
                clear spikesAll;
        end
        
        %% SPLIT DATA BACK INTO TRIALS
                
        ntrials = size(files_all,2);
        freq = []; % Initialize in case not loaded
        
        tf = strcmp(method,{'fmm','ums','ops'});
        tf = max(tf(:));
        
        if (tf) % CHANGE SO THAT ALL TRIALS ARE SAVED TO SAME OUTPUT FILE
            N = 0;
            for iTrial = 1:ntrials
                load(files_all{iProbe,iTrial});
                nspikes = size(spikes.spiketimes,2);
                assigns = assigns_all(N(iTrial)+1:N(iTrial+1));
                spikes.assigns = assigns;
                spikes   = psr_sort_clusters(spikes);
                spikes   = psr_sst_filter_rpv(spikes);
                clusters = psr_sst_clusterfeatures(spikes,freq,parameters);
                saveFile(spikes,clusters,metadata,freq,parameters,savePath); % also save LFP data
                N = N + nspikes;
            end
        else
            spikesAll    = psr_sst_filter_rpv(spikesAll,parameters); % Remove spikes based on RPVs
            % Initialize arrays
            trials       = zeros(size(spikesAll.spiketimes));
            stimVoltsAll = zeros(ntrials,1);
            onsetTimes   = zeros(ntrials,1);
            offsetTime   = 0;
            offsetSpike  = 0;
            % Store data per trial
            for iTrial = 1:ntrials
                load(files_all{iProbe,iTrial}); % Load temporary trial data
                if (iTrial == 1); stimTimesAll = cell(ntrials,length(spikes.stimtimes)); end
                stimTimesAll(iTrial,:) = spikes.stimtimes;
                stimVoltsAll(iTrial)   = metadata.stimulus;
                if (~isempty(freq))
                    freqAll(iTrial) = freq; %#ok
                end
                onsetTime   = offsetTime;
                offsetTime  = onsetTime + (size(spikes.data,2) / spikes.params.Fs);
                onsetSpike  = offsetSpike + 1;
                offsetSpike = find(spikesAll.spiketimes <= offsetTime,1,'last');
                onsetTimes(iTrial) = onsetTime;
                trials(onsetSpike:offsetSpike) = iTrial;
            end
            spikesAll.info.stimtimes  = stimTimesAll;
            spikesAll.trials          = trials;
            spikesAll.info.trialonset = onsetTimes;
            metadata.stimulus         = stimVoltsAll;
            spikesAll = psr_sst_clustermerge(spikesAll,freqAll,parameters); % Merge clusters
            saveFile(spikesAll,metadata,freqAll,parameters,savePath); % Save to MAT file 
        end
        
        %% Delete temporary spikes files
        for iTrial = 1:ntrials
            delete(files_all{iProbe,iTrial});
        end
    end
    
    delete([savePath tempFileName '.mat']);
    disp(['Spike sorting completed. MAT file(s) saved to: "' savePath '"']);
    
    
    
end

end

function saveFile(spikes,metadata,freq,parameters,savePath)

filename = [savePath ...
    'Spikes_'   metadata.subject ...
    '_'         metadata.session ...
    '_P'        num2str(metadata.probe,             '%02d') ...
    '_'         parameters.sorting.method];

if (length(metadata.stimulus) == 1)
    filename = [filename ...
        '_V'    num2str(metadata.stimulus,            '%03d')];
end

freq       = orderfields(freq);       %#ok
parameters = orderfields(parameters); %#ok
spikes     = orderfields(spikes);     %#ok

save([filename '.mat'],'spikes','metadata','freq','parameters');

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