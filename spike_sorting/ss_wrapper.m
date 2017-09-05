function ss_wrapper(subject,session,loadPath,savePath,method)

pattern = 'CH';
ext     = '.continuous';

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
    if (loadPath{iTrial}(end) ~= '\' && ~isempty(loadPath{iTrial})); loadPath{iTrial} = [loadPath{iTrial}, '\']; end
end

if nargin < 4; savePath = []; end % save in current working directory

if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

%% create metadata variable: to be expanded in future versions with data
% from the electronics notebook
metadata.tetrode = [];
metadata.subject = subject;
metadata.session = session;
metadata.dir     = cd;

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
bp_high    = 6000; % Hz
bp_low     = 1000;  % Hz
bp_order   = 10;
std_thresh = 8;
MFAtimes   = cell(numtrials,2);

pattern     = 'ADC';
fileIndices = [5 7];
nfiles      = length(fileIndices);

for iTrial = 1:numtrials
    
    if (stimuliConditions(iTrial) > 0)
    
        files = dir([loadPath{iTrial} '\*' pattern '*' ext]);
        files = char(files.name);
        files = files(fileIndices,:);
        files = strtrim(files);
        
        data_channels =  cell(nfiles,1);
        thresholds    = zeros(nfiles,1);
        
        for iFile = 1:nfiles % Filter
            
            file = [loadPath{iTrial} files(iFile,:)];
            [data, ~, info] = load_open_ephys_data(file); % data in microvolts

            Fs = info.header.sampleRate;

            % Band-pass filter

            [B,A] = butter(bp_order,[bp_low bp_high]/(Fs/2),'bandpass');
            data  = filtfilt(B,A,data);
            
            % Normalize
            
            data = data / max(data);
            
            stdev = median(abs(data)) / 0.6745; % median absolute deviation
            thresholds(iFile) = std_thresh * stdev;
            
            data_channels{iFile} = single(data); % convert to single type
            
        end
        
        threshold = max(thresholds);
        
        spikes = [];
        spikes = ss_default_params(spikes);
        spikes.params.Fs = Fs;
        [MFAtimes{iTrial,1},MFAtimes{iTrial,2}] = ss_artifact_detection_2(data_channels,threshold,spikes);
    end
    
end

%% SPIKE DETECTION / SORTING

%% user defined variables:

% Band-pass filter parameters
bp_high  = 6000; % Hz
bp_low   = 600;  % Hz
bp_order = 10;

tdata = 60; % cut data in sections of X minutes

%% per tetrode

tic

files_all = cell(1,numtrials);
ntets     = zeros(numtrials,1);
nspikes   = zeros(1,numtrials);

for iTrial = 1:numtrials
    
    data_trial  = [];
    files_trial = files_trials{iTrial};
    numfiles    = length(files_trial);
    ntets(iTrial) = numfiles / 4;
    
    for iTetrode = 1:ntets(iTrial);
                
        spikes = [];
        spikes = ss_default_params(spikes);
        
        for iElectrode = 1:4
            
            iFile = (iTetrode - 1) * 4 + iElectrode;
            
            filename = [loadPath{iTrial} files_trial{iFile}];
            filename = strtrim(filename);
            
            [data_channel, ~, info] = load_open_ephys_data(filename); % data in microvolts
            
            Fs = info.header.sampleRate;
            
            % Band-pass filter
            
            [B,A]        = butter(bp_order,[bp_low bp_high]/(Fs/2),'bandpass');
            data_channel = filtfilt(B,A,data_channel);
            data_channel = single(data_channel); % convert to single type
            
            if (iElectrode == 1) % new tetrode
                data = zeros(4,length(data_channel),'single');
                if (iTetrode == 1)
                    data_trial = zeros(4 * ntets(iTrial),length(data_channel),'single');
                end
            end
            
            data(iElectrode,:) = data_channel - mean(data_channel);
            
        end
        
        clear data_channel
        
        % Spike detection
        
        nsamples         = tdata * 60 * Fs;
        nsamples_total   = size(data,2);
        nsection         =  ceil(nsamples_total / nsamples);
        nsamples_section = floor(nsamples_total / nsection); % process data in sections
        
        iStart = 1;
        for iSection = 1:nsection
            data_section = data(:,iStart:iStart + nsamples_section - 1);
            data_section = {data_section'};
            spikes.params.Fs = Fs;  % Hz, sampling rate of spike data
            if (stimuliConditions(iTrial) > 0) % magnetic field artifact detection
                spikes = ss_artifact_detection(data_section{1},spikes,iStart);
            end
            spikes = ss_detect(data_section,spikes);
            iStart = iStart + nsamples_section;
        end
        
        filename = [savePath 'Spikes_' num2str(iTetrode) '_' num2str(iTrial) '.mat'];
        files_all{iTetrode,iTrial} = filename;
        save(filename,'spikes');
        
        data_trial(4*iTetrode-3:4*iTetrode,:) = data;
        
    end
    
    % Magnetic field artifact combination
    
    MFAtimes_1 = [];
    MFAtimes_2 = [];
    if (stimuliConditions(iTrial) > 0)
        [MFAtimes_1,MFAtimes_2] = ss_artifact_combine(files_all(:,iTrial),MFAtimes{iTrial,1},MFAtimes{iTrial,2});
    end
    
    % Artifact removal
    n = ss_artifact_removal(files_all(:,iTrial), data_trial, MFAtimes_1, MFAtimes_2);
    nspikes(1:length(n),iTrial) = n;
    
end

clear data data_section spikes

%% Clustering

%     saveFile(spikes,clusters,metadata,savePath)

ntets   = size(files_all,1);
ntrials = size(files_all,2);

for iTetrode = 1:ntets
    
    N   = sum(nspikes,2);
    itr = 1;
    spikesAll = []; % structure to contain all 'spikes' structures from all trials
    
    % Merge data across trials
    
    for iTrial = 1:ntrials
        if (strcmp(method,'kls')) % kilosort
            
        else
            load(files_all{iTetrode,iTrial});
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
    end
    
    % Unmerge data back into trials
    
    N = cumsum(nspikes,2);
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
        
        % Pre-merge save
        
        clusters = ss_clusterfeatures(spikes); 
        
        saveFile(spikes,clusters,metadata,savePath,'_1');
        
%         % Merge
%         
%         [spikes,clusters] = ss_clustermerge(spikes);
%         [spikes,clusters] = ss_spikefilter(spikes,clusters); 
%         
%         % Post-merge save
%         
%         saveFile(spikes,clusters,metadata,savePath,'_2');
        
        delete(files_all{iTetrode,iTrial}); % delete temporary spikes file
    end
    
end

toc

end

function saveFile(spikes,clusters,metadata,savePath,label) %#ok

if (nargin < 5); label = []; end

save([savePath ...
    'Spikes_'   metadata.subject ...
    '_'         metadata.session ...
    '_T'        num2str(metadata.tetrode,             '%02d') ...
    '_V'        num2str(metadata.stimulus,            '%03d') ...
    label                                                     ...
    '_'         metadata.method                               ...
    '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
    '_N'        num2str(size(spikes.waveforms,1),     '%06d') ...
    '_AG'       num2str(spikes.params.fmm_p,        '%10.2e') ...
    '.mat'],'spikes','clusters','metadata');

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