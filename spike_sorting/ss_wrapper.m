function ss_wrapper(subject,session,loadPath,savePath,method)

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

if nargin < 4;
    savePath = []; % save in current working directory
end

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

%% user defined variables:

% Band-pass filter parameters
bp_high  = 6000; % Hz
bp_low   = 600;  % Hz
bp_order = 10;

pattern = 'CH';
ext     = '.continuous';

tdata = 60; % cut data in sections of X minutes

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

files = cell(numtrials,1);

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
    files{iTrial} = files_trial;
end

%% per tetrode

tic

files_all = cell(1,numtrials);
ntets     = zeros(numtrials,1);
nspikes  = zeros(1,numtrials);

for iTrial = 1:numtrials
        
    data_trial  = [];
    files_trial = files{iTrial};
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
            spikes = ss_artifact_detection(data_section{1},spikes,iStart);
            spikes = ss_detect(data_section,spikes);
            iStart = iStart + nsamples_section;
        end
        
        filename = [savePath 'Spikes_' num2str(iTetrode) '_' num2str(iTrial) '.mat'];
        files_all{iTetrode,iTrial} = filename;
        save(filename,'spikes');
        
        data_trial(4*iTetrode-3:4*iTetrode,:) = data;
        
    end
    
    % Artifact removal
    n = ss_artifact_removal(files_all(:,iTrial), data_trial);
    nspikes(1:length(n),iTrial) = n;
    
end

clear data data_section spikes

%% Clustering

if (strcmp(method,'ums'))
    
%     spikes   = ums_clustering(spikes);
%     spikes   = ss_aggregate(spikes);
%     clusters = ums_clusterfilter(spikes); %#ok
%     
%     save([...
%         'Spikes_'   metadata.subject ...
%         '_'         metadata.session ...
%         '_T'        num2str(metadata.tetrode,             '%02d') ...
%         '_'         spikes.params.cluster_method                  ...
%         '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
%         '_N'        num2str(size(spikes.waveforms,1),     '%06d') ...
%         '_AG'       num2str(spikes.params.agg_cutoff,   '%10.2e') ...
%         '.mat'],'spikes','clusters','metadata');
    
elseif (strcmp(method,'fmm'))
    
    ntets   = size(files_all,1);
    ntrials = size(files_all,2);
    
    for iTetrode = 1:ntets
        
        N   = sum(nspikes,2);
        itr = 1;
        
        % Merge data across trials
                
        for iTrial = 1:ntrials
            load(files_all{iTetrode,iTrial});
            spikes = ss_align(spikes);
            nchan  = size(spikes.waveforms,3);
            waves  = zeros(size(spikes.waveforms,1),nchan*size(spikes.waveforms,2),'single');
            n      = size(spikes.waveforms,2);
            for ichan = 1:nchan
                waves(:,n*(ichan-1)+1:n*ichan) = spikes.waveforms(:,:,ichan);
            end
            if (iTrial == 1)
                waveforms = zeros(N(iTetrode),size(waves,2));
            end
            nspikes_trial = size(waves,1);
            waveforms(itr:itr+nspikes_trial-1,:) = waves;
            itr = itr + nspikes_trial;
        end
        
        Sorter          = FMM(waveforms',spikes.params.fmm_k);
        Sorter.align    = false;
        Sorter.FMMparam = spikes.params.fmm_p; %% Changes how aggresively to cluster, range 0-1
        Sorter.initialize;
        Sorter.runVBfit;
        
        assigns_all = getMAPassignment(Sorter);
        
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
            
            % Pre-merge save
            
            clusters = ss_clusterfeatures(spikes); %#ok
            
            save([savePath ...
                'Spikes_'   metadata.subject ...
                '_'         metadata.session ...
                '_T'        num2str(iTetrode,                     '%02d') ...
                '_V'        num2str(metadata.stimulus,            '%03d') ...
                '_1'                                                      ...
                '_'         method                                        ...
                '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
                '_N'        num2str(size(spikes.waveforms,1),     '%06d') ...
                '_AG'       num2str(spikes.params.fmm_p,        '%10.2e') ...
                '.mat'],'spikes','clusters','metadata');
            
            % Merge
            
            [spikes,clusters] = ss_clustermerge(spikes);
            [spikes,clusters] = ss_spikefilter(spikes,clusters); %#ok
            
            % Post-merge save
            
            save([savePath ...
                'Spikes_'   metadata.subject ...
                '_'         metadata.session ...
                '_T'        num2str(iTetrode,                     '%02d') ...
                '_V'        num2str(metadata.stimulus,            '%03d') ...
                '_2'                                                      ...
                '_'         method                                        ...
                '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
                '_N'        num2str(size(spikes.waveforms,1),     '%06d') ...
                '_AG'       num2str(spikes.params.fmm_p,        '%10.2e') ...
                '.mat'],'spikes','clusters','metadata');
            
            delete(files_all{iTetrode,iTrial}); % delete temporary spikes file
        end
    end
end

toc

end