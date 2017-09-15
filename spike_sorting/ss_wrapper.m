function ss_wrapper(subject,session,loadPath,savePath)

pattern = 'CH';
ext     = '.continuous';

parameters = ept_parameter_config;

% ums_wrapper ('YZ02','R170')

% metadata to include: channels
% features to include: combine different .continuous together and spike
% sort multiple sessions together; combine arbitrary sequences of
% electrodes into tetrodes; denoise (magnet artifacts);

if nargin < 3;
    loadPath = uigetdir([],'Select data folder');
end

if (~iscell(loadPath)); loadPath = {loadPath}; end

numtrials = length(loadPath);
for iTrial = 1:numtrials
    if (loadPath{iTrial}(end) ~= '\' && ~isempty(loadPath{iTrial}));
        loadPath{iTrial} = [loadPath{iTrial}, '\'];
    end
end

if nargin < 4; savePath = []; end % save in current working directory

if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

if (strcmp(parameters.cluster.method,'fmm') || strcmp(parameters.cluster.method,'ums')); SPIKE_DETECTION = 1;
else                                                                                     SPIKE_DETECTION = 0;
end

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
    filepath = fliplr(filepath);
    k = strfind(filepath,'\');
    k = k(2);
    filepath = filepath(1:k);
    filepath = fliplr(filepath);
    elems = regexp(filepath, '_(\d*)V', 'tokens', 'once');
    str   = elems{1};
    if (~isempty(str)); stimuliConditions(iTrial) = str2double(elems{1});
    else                stimuliConditions(iTrial) = 0;
    end
end

%% Find files

files_unsorted = cell(numtrials,1);
for itrial = 1:numtrials
    files_struct = dir([loadPath{itrial} '\*' pattern '*' ext]);
    if (size(files_struct,1) == 0);
        disp(['No .CONTINUOUS files in' loadPath{itrial} '... Select different path.']);
        return;
    end
    files_unsorted{itrial} = char(files_struct.name);
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

%% MAGNETIC FIELD ARTIFACT / STIMULUS ONSET DETECTION

% Band-pass filter parameters
MFAtimes   = cell(numtrials,2);

pattern     = 'ADC';
fileIndices = [5 7];
nfiles      = length(fileIndices);

for iTrial = 1:numtrials
    
    files = dir([loadPath{iTrial} '\*' pattern '*' ext]);
    files = char(files.name);
    files = files(fileIndices,:);
    files = strtrim(files);
    
    data_channels =  cell(nfiles,1);
    thresholds    = zeros(nfiles,1);
    Fs_array      = zeros(nfiles,1);
    
    for iFile = 1:nfiles % Filter
        
        file = [loadPath{iTrial} files(iFile,:)];
        [data_tetrode_sst, ~, info] = load_open_ephys_data(file); % data in microvolts
        
        Fs = info.header.sampleRate;
        Fs_array(iFile) = Fs;
        
        % Band-pass filter
        
        [B,A] = butter(parameters.mfa.bp_order,[parameters.mfa.bp_low parameters.mfa.bp_high]/(Fs/2),'bandpass');
        data_tetrode_sst  = filtfilt(B,A,data_tetrode_sst);
        
        % Normalize
        
        data_tetrode_sst = data_tetrode_sst / max(data_tetrode_sst);
        
        stdev = median(abs(data_tetrode_sst)) / 0.6745; % median absolute deviation
        thresholds(iFile) = parameters.mfa.thresh * stdev;
        
        data_channels{iFile} = single(data_tetrode_sst); % convert to single type
        
    end
    
    Fs        = mode(Fs_array);
    threshold = max(thresholds);
    
    if (stimuliConditions(iTrial) > 0)
        [MFAtimes{iTrial,1},MFAtimes{iTrial,2}] = ss_mfa_detection(data_channels,threshold,parameters,Fs);
    else % randomly select MFA times in order to get non-stimulus control data
        dur = length(data_channels{1}) / Fs; 
        MFAtimes{iTrial,1} = sort(dur * rand(1,parameters.mfa.control));
        MFAtimes{iTrial,2} = sort(dur * rand(1,parameters.mfa.control));
    end
end

%% SPIKE DETECTION / SORTING

%% user defined variables:

%% per tetrode

tic

files_all =  cell(1,numtrials);
ntets     = zeros(numtrials,1);
nlength   = zeros(1,numtrials);

Fr = 4 * parameters.lfp.bp_high; % Resample at 4 times the highest rate

for iTrial = 1:numtrials
    
    data_trial_sst = [];
    files_trial    = files_trials{iTrial};
    numfiles       = length(files_trial);
    ntets(iTrial)  = numfiles / 4;
    
    for iTetrode = 1:ntets(iTrial);
        
        for iElectrode = 1:parameters.nelectrodes
            
            % Load Open-Ephys data
            
            iFile = (iTetrode - 1) * 4 + iElectrode;
            
            filename = [loadPath{iTrial} files_trial{iFile}];
            filename = strtrim(filename);
            
            [data_channel_raw, timestamps, info] = load_open_ephys_data(filename); % data in microvolts
            
            Fs = info.header.sampleRate;
            
            % Create arrays
            
            if (iElectrode == 1);
                Fs_array         = zeros(parameters.nelectrodes,1);
                data_tetrode_sst = zeros(parameters.nelectrodes,length(data_channel_raw),'single');
                if (iTetrode == 1);
                    data_trial_sst = zeros(4 * ntets(iTrial),length(data_channel_raw),'single');
                end
            end
            
            Fs_array(iElectrode) = Fs; % Sampling frequencies should be equal
            
            %% BAND-PASS FILTERING SPIKE SORTING
            
            [B,A]        = butter(parameters.spikes.bp_order,[parameters.spikes.bp_low parameters.spikes.bp_high]/(Fs/2),'bandpass');
            data_channel = filtfilt(B,A,data_channel_raw);
            data_channel = single(data_channel); % convert to single type
            
            data_tetrode_sst(iElectrode,:) = data_channel - mean(data_channel);
            
            %% LOCAL FIELD POTENTIAL
            
            % Resample raw signal
            
            timestamps        = timestamps - timestamps(1);
            [data,timestamps] = resample(data_channel_raw,timestamps,Fr); % resample
            
            if (iElectrode == 1);
                data_tetrode_lfp = zeros(parameters.nelectrodes,length(data),'single');
                time_tetrode_lfp = zeros(parameters.nelectrodes,length(data),'single');
            end
            
            time_tetrode_lfp(iElectrode,:) = timestamps;
            data_tetrode_lfp(iElectrode,:) = data;
        end
        
        clear data_channel
        
        flag = sum(Fs_array ~= Fs_array(1)) ~= 0;
        if (flag); disp('Sampling frequency mismatch'); continue;
        else Fs = Fs_array(1);
        end
        
        % SPIKE DETECTION
        
        nsamples         = 60 * parameters.spikes.tsection * Fs; % cut data in sections of X minutes
        nsamples_total   = size(data_tetrode_sst,2);
        nsection         =  ceil(nsamples_total / nsamples);
        nsamples_section = floor(nsamples_total / nsection); % process data in sections
        
        spikes = [];
        
        if (SPIKE_DETECTION)
            spikes = ss_default_params(spikes);
            iStart = 1;
            for iSection = 1:nsection
                data_section = data_tetrode_sst(:,iStart:iStart + nsamples_section - 1);
                data_section = {data_section'};
                spikes.params.Fs = Fs;  % Hz, sampling rate of spike data
                %             if (stimuliConditions(iTrial) > 0) % magnetic field artifact detection
                %                 spikes = ss_mfa_detection_raw(data_section{1},spikes,iStart);
                %             end
                spikes = ss_detect(data_section,spikes);
                iStart = iStart + nsamples_section;
            end
        else % save raw data
            spikes.data      = data_tetrode_sst;
            spikes.params.Fs = Fs;
            spikes.nlength   = size(data_tetrode_sst,2);
        end
        
        data_trial_sst(4*iTetrode-3:4*iTetrode,:) = data_tetrode_sst;
        
        % LOCAL FIELD POTENTIAL
        
        % Maybe resample first
        
        timestamps = mean(time_tetrode_lfp,1);
        data       = ept_convert2fieldtrip(data_tetrode_lfp,timestamps,MFAtimes{iTrial,1},parameters,Fr);
        data       = ept_preprocessing(data,parameters); % FT preprocessing
        freq       = 0;
        if (~isempty(data)); freq = ept_timefreq_analysis(data,parameters); end
        
        % Save
        
        filename = [savePath 'Spikes_' num2str(iTetrode) '_' num2str(iTrial) '.mat'];
        files_all{iTetrode,iTrial} = filename;
        save(filename,'spikes','freq','parameters');
    end
    
    %     % Magnetic field artifact combination
    %
    %     MFAtimes_1 = [];
    %     MFAtimes_2 = [];
    %     if (stimuliConditions(iTrial) > 0)
    %         [MFAtimes_1,MFAtimes_2] = ss_mfa_combine(files_all(:,iTrial),MFAtimes{iTrial,1},MFAtimes{iTrial,2});
    %     end
    
    if (SPIKE_DETECTION)
        MFAtimes_1 = MFAtimes{iTrial,1}; % without raw + control MFA combination
        MFAtimes_2 = MFAtimes{iTrial,2};
        
        % Artifact removal
        n = ss_artifact_removal(files_all(:,iTrial), data_trial_sst, MFAtimes_1, MFAtimes_2);
        nlength(1:length(n),iTrial) = n;
    end
end

%% CLUSTERING
keepvars = {'parameters','files_all','nlength','stimuliConditions','metadata','savePath','SPIKE_DETECTION'};
clearvars('-except', keepvars{:});
ntets   = size(files_all,1);
ntrials = size(files_all,2);
method  = parameters.cluster.method;

for iTetrode = 1:ntets % do clustering per tetrode across all stimulus conditions
    
    N   = sum(nlength,2);
    itr = 1;
    spikesAll = []; % structure to contain all 'spikes' structures from all trials
    
    data_tetrode      = []; %
    data_tetrode.data = [];
    
    % Merge data across trials
    
    for iTrial = 1:ntrials
        load(files_all{iTetrode,iTrial});
        if (~SPIKE_DETECTION) % method requires raw data
            if (strcmp(method,'cbp'))
                data_tetrode.data = [data_tetrode.data,spikes.data];
                data_tetrode.dt   = 1 / spikes.params.Fs;
            elseif (strcmp(method,'kls'))
                
            end
        else % method uses detected spike waveforms
            spikes = ss_align(spikes);
            if (strcmp(method,'fmm')) % dictionary learning
                % Convert waveforms for dictionary learning
                nchan  = size(spikes.waveforms,3);
                waves  = zeros(size(spikes.waveforms,1),nchan*size(spikes.waveforms,2),'single');
                n      = size(spikes.waveforms,2);
                for ichan = 1:nchan
                    waves(:,n*(ichan-1)+1:n*ichan) = spikes.waveforms(:,:,ichan);
                end
                if (iTrial == 1); waveforms = zeros(N(iTetrode),size(waves,2)); end
                nspikes_trial = size(waves,1);
                waveforms(itr:itr+nspikes_trial-1,:) = waves;
                itr = itr + nspikes_trial;
            elseif (strcmp(method,'ums')) % UltraMegaSort
                spikesAll = appendSpikes(spikesAll,spikes);
            end
        end
    end
    
    if (strcmp(method,'ums')) % UltraMegaSort
        spikesAll = meanSpikes  (spikesAll);
        spikesAll = ss_kmeans   (spikesAll);
        spikesAll = ss_energy   (spikesAll);
        spikesAll = ss_aggregate(spikesAll);
        assigns_all = spikesAll.assigns;
        clear spikesAll;
    elseif (strcmp(method,'fmm')) % dictionary learning
        assigns_all = ss_dictionary_learning(spikes,waveforms);
    elseif (strcmp(method,'cbp')) % bayesian persuit
        data_tetrode.nchan    = size(data_tetrode.data,1);
        data_tetrode.nsamples = size(data_tetrode.data,2);
        data_tetrode.polarity = 'min';
        spike_sorting_cbp(data_tetrode);
    elseif (strcmp(method,'kls')) % kilosort
        
    end
    
    % Unmerge data back into trials
    
    N = cumsum(nlength,2);
    N = [0,N(iTetrode,:)];
    
    for iTrial = 1:size(files_all,2)
        
        load(files_all{iTetrode,iTrial});
        assigns = assigns_all(N(iTrial)+1:N(iTrial+1));
        
        spikes.assigns             = assigns;
        spikes.info.kmeans.assigns = assigns;
        numclusts                  = max(spikes.info.kmeans.assigns);
        cmap                       = jetm(numclusts);
        spikes.info.kmeans.colors  = cmap(randperm(numclusts),:);
        
        metadata.stimulus = stimuliConditions(iTrial);
        metadata.tetrode  = iTetrode;
        metadata.method   = method;
        
        spikes = ss_spikefilter_rpv(spikes);
        clusters = ss_clusterfeatures(spikes);
        saveFile(spikes,clusters,metadata,freq,parameters,savePath); % also save LFP data
        
        %         % Merge
        %         [spikes,clusters] = ss_clustermerge(spikes);
        %         % Post-merge save
        %         saveFile(spikes,clusters,metadata,savePath,'_2');
        
        delete(files_all{iTetrode,iTrial}); % delete temporary spikes file
    end
    
end

toc

end

function saveFile(spikes,clusters,metadata,freq,parameters,savePath,label) %#ok

if (nargin < 7); label = []; end

if (strcmp(parameters.cluster.method,'fmm')); agg_parameter = spikes.params.sorting.fmm.p;
else                                          agg_parameter = [];
end

save([savePath ...
    'Spikes_'   metadata.subject ...
    '_'         metadata.session ...
    '_T'        num2str(metadata.tetrode,             '%02d')   ...
    '_V'        num2str(metadata.stimulus,            '%03d')   ...
    label                                                       ...
    '_'         metadata.method                                 ...
    '_ST'       num2str(spikes.params.detect.thresh,  '%04.1f') ...
    '_N'        num2str(size(spikes.waveforms,1),     '%06d')   ...
    '_AG'       num2str(agg_parameter,                '%10.2e') ...
    '.mat'],'spikes','clusters','metadata','freq','parameters');

end

function spikesMain = appendSpikes(spikesMain,spikes)

spikesMain.waveforms       = [spikesMain.waveforms       ; spikes.waveforms];
spikesMain.spiketimes      = [spikesMain.spiketimes      ; spikes.spiketimes];
spikesMain.trials          = [spikesMain.trials          ; spikes.trials];
spikesMain.nspikes         = [spikesMain.nspikes         ; spikes.nspikes];
spikesMain.unwrapped_times = [spikesMain.unwrapped_times ; spikes.unwrapped_times];
spikesMain.assigns         = [spikesMain.assigns         ; spikes.assigns];

spikesMain.info.detect.stds   = [spikesMain.info.detect.stds   ; spikes.info.detect.stds];
spikesMain.info.detect.thresh = [spikesMain.info.detect.thresh ; spikes.info.detect.thresh];
end

function spikes = meanSpikes(spikes)

spikes.info.detect.stds   = mean(spikes.info.detect.stds);
spikes.info.detect.thresh = mean(spikes.info.detect.thresh);

end