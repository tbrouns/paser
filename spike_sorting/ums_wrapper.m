function ums_wrapper(subject,session,loadPath,savePath)

% ums_wrapper ('YZ02','R170')

% metadata to include: channels
% features to include: combine different .continuous together and spike
% sort multiple sessions together; combine arbitrary sequences of
% electrodes into tetrodes; denoise (magnet artifacts);

if nargin < 3;
    loadPath = uigetdir([],'Select data folder');
end

if nargin < 4;
    savePath = []; % save in current working directory
end

if (loadPath(end) ~= '\' && ~isempty(loadPath)); loadPath = [loadPath, '\']; end
if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

%% create metadata variable: to be expanded in future versions with data
% from the electronics notebook
metadata.tetrode = [];
metadata.subject = subject;
metadata.session = session;
metadata.dir     = cd;

clear subject session

%% user defined variables:

% Sets termination criterion for cluster aggregation. Higher values allow
% less aggregation. Lower values allow more aggregation. Aggregation is
% stopped when overlap density is less than the given % of main cluster
% density

% Band-pass filter parameters
bp_high  = 6000; % Hz
bp_low   = 600;  % Hz
bp_order = 10;

pattern = 'CH';
ext     = '.continuous';

tdata = 60; % cut data in sections of X minutes

%% Find files

files_unsorted = dir([loadPath '\*' pattern '*' ext]);
if (size(files_unsorted,1) == 0);
    disp('No .CONTINUOUS files in folder. Select different path.');
    return;
end
files_unsorted = char(files_unsorted.name);

%% sort files

numfiles = length(files_unsorted(:,1));

files = cell(numfiles,1);
for iFile = 1:numfiles
    filename   = files_unsorted(iFile,:);
    filename   = strtrim(filename);
    k          = strfind(filename,pattern) + length(pattern);
    [~,name,~] = fileparts(filename);
    id         = str2double(name(k:end));
    files{id}  = filename;
end

files = files(~cellfun('isempty',files)); % remove empty cells

%% per tetrode

numtets = numfiles / 4;

for iTetrode = 1:numtets;
    tic
    metadata.tetrode = iTetrode;
    for iElectrode = 1:4; % load all tetrode data
        
        iFile = (iTetrode - 1) * 4 + iElectrode;
        
        filename = [loadPath files{iFile}];
        filename = strtrim(filename);
        
        [data_channel, ~, info] = load_open_ephys_data(filename); % data in microvolts
        
        Fs  = info.header.sampleRate;
        
        % Band-pass filter
        
        [B,A]        = butter(bp_order,[bp_low bp_high]/(Fs/2),'bandpass');
        data_channel = filtfilt(B,A,data_channel);
        data_channel = single(data_channel); % convert to single type
        
        if (iElectrode == 1) % new tetrode
            labels  = cell(4,1);
            data    = zeros(4,length(data_channel),'single');
        end
        
        labels{iElectrode} = num2str(iFile);
        data(iElectrode,:) = data_channel - mean(data_channel);
        
    end
    
    clear data_channel
    
    spikes = ss_default_params(Fs); % set parameters for spike detection
        
    % UMS spike detection
    
    nsamples_section = tdata * 60 * Fs;
    nsamples_total   = size(data,2);
    nsection         =  ceil(nsamples_total / nsamples_section);
    nsamples_section = floor(nsamples_total / nsection); % process data in sections
    iStart           = 1;
        
    for iSection = 1:nsection
        data_section = data(:,iStart:iStart + nsamples_section - 1);
        data_section = {data_section'};
        if (spikes.params.artifact_removal); 
            [data_section,spikes] = ums_artifact_removal(data_section,spikes); 
        end
        spikes       = ss_detect(data_section,spikes);
        iStart       = iStart + nsamples_section;
    end
    
    %clear data data_section
    
    % Cluster detected spikes
    
    if (strcmp(spikes.params.cluster_method,'ums'))
        
        spikes   = ums_clustering(spikes);
        spikes   = ss_aggregate(spikes);
        clusters = ums_clusterfilter(spikes); %#ok
        
        save([... 
        'Spikes_'   metadata.subject ...
        '_'         metadata.session ...
        '_T'        num2str(metadata.tetrode,             '%02d') ...
        '_'         spikes.params.cluster_method                  ...
        '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
        '_N'        num2str(size(spikes.waveforms,1),     '%06d') ...
        '_AG'       num2str(spikes.params.agg_cutoff,   '%10.2e') ...
        '.mat'],'spikes','clusters','metadata');
        
    elseif (strcmp(spikes.params.cluster_method,'fmm'))
        
        spikes    = ss_align(spikes);
        nchan     = size(spikes.waveforms,3);
        waveforms = zeros(size(spikes.waveforms,1),nchan*size(spikes.waveforms,2));
        n         = size(spikes.waveforms,2);
        for ichan = 1:nchan
            waveforms(:,n*(ichan-1)+1:n*ichan) = spikes.waveforms(:,:,ichan);
        end
        
        Sorter          = FMM(waveforms',spikes.params.fmm_k);
        Sorter.align    = false;
        Sorter.FMMparam = spikes.params.fmm_p; %% Changes how aggresively to cluster, range 0-1
        Sorter.initialize;
        Sorter.runVBfit;
        
        assigns                    = getMAPassignment(Sorter);
        spikes.assigns             = assigns;
        spikes.info.kmeans.assigns = assigns;
        numclusts                  = max(spikes.info.kmeans.assigns);
        cmap                       = jetm(numclusts);
        spikes.info.kmeans.colors  = cmap(randperm(numclusts),:);
        
        clusters = ums_clusterfilter(spikes); %#ok
        
        save([savePath ... 
        'Spikes_'   metadata.subject ...
        '_'         metadata.session ...
        '_T'        num2str(metadata.tetrode,             '%02d') ...
        '_'         spikes.params.cluster_method                  ...
        '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
        '_N'        num2str(size(spikes.waveforms,1),     '%06d') ...
        '_AG'       num2str(spikes.params.fmm_p,        '%10.2e') ...
        '.mat'],'spikes','clusters','metadata');
    end
        
    toc
    
end

end