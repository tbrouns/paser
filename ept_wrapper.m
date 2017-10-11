function ept_wrapper(subject,session,loadPath,savePath)

% fig = figure; set(gcf,'position',get(0,'screensize')); % TEMP

parameters   = ept_parameter_config;
nelectrodes  = parameters.general.nelectrodes;
pattern      = parameters.general.pattern;
ext          = parameters.general.ext;
method       = lower(parameters.sorting.method);
tempFileName = 'Temp_vars';

% ums_wrapper ('YZ02','R170')

% metadata to include: channels
% features to include: combine different .continuous together and spike
% sort multiple sessions together; combine arbitrary sequences of
% electrodes into tetrodes; denoise (magnet artifacts);

if nargin < 3; loadPath = uigetdir([],'Select data folder'); end
if nargin < 4; savePath = []; end % save in current working directory

if (~iscell(loadPath)); loadPath = {loadPath}; end

numtrials = length(loadPath);
for iTrial = 1:numtrials
    if (loadPath{iTrial}(end) ~= '\' && ~isempty(loadPath{iTrial}))
        loadPath{iTrial} = [loadPath{iTrial}, '\'];
    end
end

if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

% Check if spike sorting method requires spike detection
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
files_temp = dir([savePath 'Temp_Tetrode_*_Trial_*.mat']);
files_temp = char(files_temp.name);
nfiles     = size(files_temp,1);
for iFile = 1:nfiles
    filename = [savePath files_temp(iFile,:)];
    tetrode  = extractStringFromPath(filename,'Tetrode_(\d*)_');
    trial    = extractStringFromPath(filename,'Trial_(\d*).');
    files_all{tetrode,trial} = filename;
end
file_temp = dir([savePath tempFileName '.mat']);
file_temp = char(file_temp.name);

if (isempty(files_all) || isempty(file_temp))
    
    %% create metadata variable: to be expanded in future versions with data
    % from the electronics notebook
    
    metadata = [];
    metadata.subject  = subject;
    metadata.session  = session;
    metadata.dir      = cd;
    clear subject session
    
    % Get stimulus conditions
    
    stimuliConditions = zeros(numtrials,1);
    
    for iTrial = 1:numtrials
        filepath = loadPath{iTrial};
        stimuliConditions(iTrial) = extractStringFromPath(filepath,'_(\d*)V');
    end
    
    %% Find files
    
    files_unsorted = cell(numtrials,1);
    for iTrial = 1:numtrials
        files_struct = dir([loadPath{iTrial} '\*' pattern '*' ext]);
        if (size(files_struct,1) == 0); return; end
        files_unsorted{iTrial} = char(files_struct.name);
    end
    
    %% sort files
    
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
    
    %% STIMULUS ONSET DETECTION (MAGNETIC FIELD ARTIFACTS)
    
    MFAtimes = cell(numtrials,2);
    
    if (parameters.process.mfa)
        
        % Band-pass filter parameters
        
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
            
            data_channels =  cell(nfiles,1);
            thresholds    = zeros(nfiles,1);
            Fs_array      = zeros(nfiles,1);
            
            for iFile = 1:nfiles % Filter
                
                file = [loadPath{iTrial} files(iFile,:)];
                [data_ADC, ~, info] = load_open_ephys_data(file); % data in microvolts
                
                Fs = info.header.sampleRate;
                Fs_array(iFile) = Fs;
                
                % Band-pass filter
                
                [B,A] = butter(parameters.mfa.bp_order,[parameters.mfa.bp_low parameters.mfa.bp_high]/(Fs/2),'bandpass');
                data_ADC = filtfilt(B,A,data_ADC); % Zero-phase digital filtering
                
                % Normalize
                
                data_ADC = data_ADC / max(data_ADC);
                
                stdev = median(abs(data_ADC)) / 0.6745; % median absolute deviation
                thresholds(iFile) = parameters.mfa.thresh * stdev;
                
                data_channels{iFile} = single(data_ADC); % convert to single type
                
            end
            
            Fs        = mode(Fs_array);
            threshold = max(thresholds);
            
            if (stimuliConditions(iTrial) > 0)
                [MFAtimes{iTrial,1},MFAtimes{iTrial,2}] = ept_stim_detection(data_channels,threshold,parameters,Fs);
            else % randomly select MFA times in order to get non-stimulus control data
                dur = length(data_channels{1}) / Fs;
                MFAtimes{iTrial,1} = sort(dur * rand(1,parameters.mfa.control));
                MFAtimes{iTrial,2} = sort(dur * rand(1,parameters.mfa.control));
            end
        end
        
    end
    
    %% SPIKE + LFP DETECTION
    
    %% Per tetrode
    
    ntets   = zeros(numtrials,1);
    
    Fr = parameters.lfp.rsfactor * parameters.lfp.bp_high; % Resample at 4 times the highest rate
    
    for iTrial = 1:numtrials
        
        % Check for existing files
        
        if (iTrial <= size(files_all,2))
            files_trial = cell2mat(files_all(:,iTrial));
            numfiles    = size(files_trial,1);
            ntets       = size(files_trials{iTrial},1) / parameters.general.nelectrodes;
            if (numfiles == ntets); continue; end
        end
        
        % Initialize
        
        files_trial   = files_trials{iTrial};
        numfiles      = length(files_trial);
        ntets(iTrial) = numfiles / parameters.general.nelectrodes;
        MFAtimesTrial = MFAtimes(iTrial,:);
        stim          = stimuliConditions(iTrial);
        
        for iTetrode = 1:ntets(iTrial)
            
            for iElectrode = 1:nelectrodes
                
                % Load Open-Ephys data
                
                iFile = (iTetrode - 1) * nelectrodes + iElectrode;
                
                filename = [loadPath{iTrial} files_trial{iFile}];
                filename = strtrim(filename);
                
                [data_channel_raw, timestamps, info] = load_open_ephys_data(filename); % data in microvolts
                
                % Store sampling frequency
                Fs = info.header.sampleRate;
                if (iElectrode == 1); Fs_array = zeros(nelectrodes,1); end
                Fs_array(iElectrode) = Fs; % Sampling frequencies should be equal
                
                data_channel_raw = ept_artifact_fft(data_channel_raw,parameters,Fs);
                
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
%                 export_fig([savePath 'FFT_' num2str(iTrial) '_' num2str(iTetrode) '_' num2str(iElectrode)]);
                %%%%
                
                if (parameters.process.spikes)
                    
                    %% Initialize arrays
                    if (iElectrode == 1)
                        data_tetrode_sst = zeros(nelectrodes,length(data_channel_raw),'single');
                        if (iTetrode == 1 && parameters.spikes.artifacts_removal)
                            data_trial_sst = zeros(ntets(iTrial)*nelectrodes,length(data_channel_raw),'single');
                        end
                    end
                    
                    %% BAND-PASS FILTERING SPIKE SORTING
                    
                    [B,A]        = butter(parameters.spikes.bp_order,[parameters.spikes.bp_low parameters.spikes.bp_high]/(Fs/2),'bandpass');
                    data_channel = filtfilt(B,A,data_channel_raw); % Zero-phase digital filtering
                    data_channel = single(data_channel); % convert to single type
                    
                    data_tetrode_sst(iElectrode,:) = detrend(data_channel, 'linear');
                    
                end
                
                %% LOCAL FIELD POTENTIAL
                
                % Resample raw signal
                
                timestamps        = timestamps - timestamps(1);
                [data,timestamps] = resample(data_channel_raw,timestamps,Fr); % resample
                
                if (iElectrode == 1)
                    data_tetrode_lfp = zeros(nelectrodes,length(data),'single');
                    time_tetrode_lfp = zeros(nelectrodes,length(data),'single');
                end
                
                time_tetrode_lfp(iElectrode,:) = timestamps;
                data_tetrode_lfp(iElectrode,:) = data;
            end
            
            %%%% TEMP
            %             figure(fig);
            %             plot((1:size(data_tetrode_sst,2))/Fs,data_tetrode_sst);
            %             xlabel('Time [s]');
            %             ylabel('Voltage');
            %             savefig([savePath 'Signal_' num2str(iTrial) '_' num2str(iTetrode) '.fig']);
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
                
                nsamples         = 60 * parameters.spikes.tsection * Fs; % cut data in sections
                nsamples_total   = size(data_tetrode_sst,2);
                nsection         =  ceil(nsamples_total / nsamples);
                nsamples_section = floor(nsamples_total / nsection); % process data in sections
                
                if (SPIKE_DETECTION)
                    iStart = 1;
                    for iSection = 1:nsection
                        data_section = data_tetrode_sst(:,iStart:iStart + nsamples_section - 1);
                        data_section = {data_section'};
                        if (stim > 0 && parameters.spikes.artifacts_combine) % magnetic field artifact detection
                            spikes = ss_mfa_detection_raw(data_section{1},spikes,iStart);
                        end
                        spikes = ept_sst_detect(data_section,spikes);
                        iStart = iStart + nsamples_section;
                    end
                end
                
                if (SAVE_RAW_DATA)
                    spikes.data = data_tetrode_sst;
                end
                
                if (parameters.spikes.artifacts_removal)
                    data_trial_sst(nelectrodes*(iTetrode-1)+1:nelectrodes*iTetrode,:) = data_tetrode_sst;
                end
            end
            
            spikes.stimtimes = MFAtimesTrial;
            
            % LOCAL FIELD POTENTIAL
            
            freq = [];
            if (parameters.process.lfp && ~isempty(spikes.stimtimes{1})) % TEMP: also demand that stimtimes are avaiable (change for active data)
                timestamps = mean(time_tetrode_lfp,1);
                [data,artifacts] = ept_artifact_removal_lfp(data_tetrode_lfp,Fr,parameters);
                data = ept_convert2fieldtrip(data,timestamps,MFAtimes{iTrial,1},parameters,Fr);
                data = ept_preprocessing(data,parameters); % FT preprocessing
                if (~isempty(data)); freq = ept_timefreq_analysis(data,parameters,stim); end
                freq.artifacts = artifacts;
            end
            
            % Save
            
            metadata.stimulus = stimuliConditions(iTrial);
            metadata.tetrode  = iTetrode;
            
            if (parameters.process.spikes)
                % Save temporary MAT file
                filename = [savePath 'Temp_Tetrode_' num2str(iTetrode) '_Trial_' num2str(iTrial) '.mat'];
                files_all{iTetrode,iTrial} = filename;
                save(filename,'spikes','metadata','freq','parameters');
            else
                % Save LFP output
                saveFile(spikes,[],metadata,freq,parameters,savePath);
            end
        end
        
        %% Filter spiking data across all tetrodes
                        
        if (parameters.process.spikes)
            
            filesTrial = files_all(:,iTrial);
            
            if (stim > 0 && parameters.spikes.artifacts_combine)
                MFAtimesTrial = ss_mfa_combine(filesTrial,MFAtimesTrial); % Magnetic field artifact combination
                ept_mfa_save(filesTrial,MFAtimesTrial)
            end
            
            %% Global median subtraction
            if (parameters.spikes.artifacts_subtract)
                ept_sst_artifact_subtraction(filesTrial, data_trial_sst);
            end
            
            %% Artifact removal
            if (SPIKE_DETECTION && parameters.spikes.artifacts_removal)
                ept_sst_artifact_removal(filesTrial, data_trial_sst, MFAtimesTrial);
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
    ntets   = size(files_all,1);
    ntrials = size(files_all,2);
    
    for iTetrode = 1:ntets % Do clustering per tetrode across all stimulus conditions
        
        % Check data
        
        MISSING_FILE = false;
        for iTrial = 1:ntrials
            if (~exist(files_all{iTetrode,iTrial},'file')); MISSING_FILE = true; break; end
        end
        if (MISSING_FILE) % If not all files are present
            for iTrial = 1:ntrials
                filename = files_all{iTetrode,iTrial};
                if (exist(filename,'file')); delete(filename); end % Delete temporary spikes file
            end
            disp(['Skipping tetrode ' num2str(iTetrode) '...']);
            continue; % Move on to next tetrode
        end
        
        disp(['Spike sorting tetrode ' num2str(iTetrode) '...']);
        
        % Initialize
        
        switch method
            case 'cbp'
                data_tetrode      = []; %
                data_tetrode.data = [];
            case {'fmm','ops','ums'}
                spikesAll = []; % structure to contain all 'spikes' structures from all trials
            case 'kst'
                data_tetrode = [];
        end
        
        %% MERGE DATA ACROSS TRIALS
        
        for iTrial = 1:ntrials
            
            load(files_all{iTetrode,iTrial});
            
            switch method
                case 'cbp'
                    data_tetrode.data = [data_tetrode.data,spikes.data];
                    data_tetrode.dt   = 1 / spikes.params.Fs;
                case 'fmm'
                    spikesAll = ept_sst_spike_append(spikesAll,spikes);
                case 'kst'
                    data_tetrode = [data_tetrode,spikes.data]; %#ok
                    Fs = spikes.params.Fs;
                case 'ops'
                    spikesAll = ept_sst_spike_append(spikesAll,spikes);
                case 'ums'
                    spikesAll = ept_sst_spike_append(spikesAll,spikes);
            end
        end
        
        data_tetrode = real(data_tetrode); % TEMP
        
        %% SPIKE SORT
        
        switch method
            case 'cbp' % WORK IN PROGRESS
                assigns_all = ept_sst_sorting_CBP(data_tetrode,parameters);
                clear data_tetrode;
            case 'fmm'
                assigns_all = ept_sst_sorting_FMM(spikesAll,parameters);
            case 'kst'
                parameters.Fs = Fs;
                rez = ept_sst_sorting_KST(data_tetrode,parameters,savePath);
                spikesAll = ept_kst_convert2spikes(rez,data_tetrode,parameters);
                spikesAll.params.Fs = rez.ops.fs;
                delete(rez.ops.fbinary);
                clear data_tetrode;
            case 'ops' % WORK IN PROGRESS
                assigns_all = ept_sst_sorting_OPS(spikesAll,parameters);
                clear spikesAll;
            case 'ums'
                assigns_all = ept_sst_sorting_UMS(spikesAll,parameters);
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
                load(files_all{iTetrode,iTrial});
                nspikes = size(spikes.spiketimes,2);
                assigns = assigns_all(N(iTrial)+1:N(iTrial+1));
                spikes.assigns = assigns;
                spikes   = ept_sort_clusters(spikes);
                spikes   = ept_sst_filter_rpv(spikes);
                clusters = ept_sst_clusterfeatures(spikes,freq,parameters);
                saveFile(spikes,clusters,metadata,freq,parameters,savePath); % also save LFP data
                N = N + nspikes;
            end
        else
            parameters = ept_parameter_config(); % TEMP
            spikesAll  = ept_sst_filter_rpv(spikesAll,parameters);
            trials       = zeros(size(spikesAll.spiketimes));
            stimVoltsAll = zeros(ntrials,1);
            onsetTimes  = zeros(ntrials,1);
            offsetTime  = 0;
            offsetSpike = 0;
            for iTrial = 1:ntrials
                load(files_all{iTetrode,iTrial});
                if (iTrial == 1); stimTimesAll = cell(ntrials,length(spikes.stimtimes)); end
                stimTimesAll(iTrial,:) = spikes.stimtimes;
                stimVoltsAll(iTrial)   = metadata.stimulus;
                freqAll(iTrial)        = freq; %#ok
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
            spikesAll = ept_sst_clustermerge(spikesAll,freqAll,parameters);
            saveFile(spikesAll,metadata,freqAll,parameters,savePath); % also save LFP data
        end
        
        %% Delete temporary spikes files
        for iTrial = 1:ntrials
            delete(files_all{iTetrode,iTrial});
        end
    end
    
    delete([savePath tempFileName '.mat']);
    disp(['Spike sorting completed. Data is saved to: "' savePath '"']);
    
end

end

function saveFile(spikes,metadata,freq,parameters,savePath) %#ok

filename = [savePath ...
    'Spikes_'   metadata.subject ...
    '_'         metadata.session ...
    '_T'        num2str(metadata.tetrode,             '%02d') ...
    '_'         parameters.sorting.method];

if (length(metadata.stimulus) == 1)
    filename = [filename ...
        '_V'    num2str(metadata.stimulus,            '%03d')];
end

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