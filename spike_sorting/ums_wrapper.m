function ums_wrapper(subject,session,fpath)

% ums_wrapper ('YZ02','R170')

% metadata to include: channels
% features to include: combine different .continuous together and spike
% sort multiple sessions together; combine arbitrary sequences of
% electrodes into tetrodes; denoise (magnet artifacts); 


if nargin < 3;
    fpath = [];
end

%% create metadata variable: to be expanded in future versions with data
% from the electronics notebook
metadata.tetrode = [];
metadata.subject = subject;
metadata.session = session;
metadata.dir = cd;

clear subject session
%% user defined variables:

% Sets termination criterion for cluster aggregation. Higher values allow
% less aggregation. Lower values allow more aggregation.
agg_cutoff = [0.00001 0.0001 0.001]; 

refractory_period = 1.5; % ms
threshold         = 4.0; % stds

detect_method = 'mad'; % maximum absolute deviation

% Band-pass filter parameters
bp_high     = 6000;
bp_low      = 600;
bp_order    = 10;

pattern = 'CH';
ext     = '.continuous';

tdata = 10; % cut data in sections of X minutes

%% Find files

files_unsorted    = dir([fpath '*' pattern '*' ext]);
if (size(files_unsorted,1) == 0);
    disp('No .CONTINUOUS files in folder. Select different path.');
    return;
end
files_unsorted    = char(files_unsorted.name);

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
    for iElectrode = 1:4; % load all tetrode data
        
        iFile = (iTetrode - 1) * 4 + iElectrode;
        
        filename = [fpath files{iFile}];
        filename = strtrim(filename);
        
        [data_channel, ~, info] = load_open_ephys_data(filename); % data in microvolts
        
        Fs  = info.header.sampleRate;
        
        if (iElectrode == 1) % new tetrode
            labels  = cell(4,1);
            data    = zeros(4,length(data_channel));
        end
        
        % Band-pass filter
        
        [B,A]           = butter(bp_order,[bp_low bp_high]/(Fs/2),'bandpass');
        data_channel    = filtfilt(B,A,data_channel);
        
        labels{iElectrode} = num2str(iFile);
        data(iElectrode,:) = data_channel;
        
    end
    
    clear data_channel
    
    for AGcutoff = 1:length(agg_cutoff)
        
        % set parameters for spike detection
        
        spikes = ss_default_params(Fs);
        spikes.params.agg_cutoff        = agg_cutoff(AGcutoff);
        spikes.params.detect_method     = detect_method;
        spikes.params.refractory_period = refractory_period;
        spikes.params.thresh            = threshold;
        
        % UMS spike detection
        
        nsamples_section = tdata * 60 * Fs;
        nsamples_total   = size(data,2);
        nsection         = floor(nsamples_total / nsamples_section);
        nsamples_section = floor(nsamples_total / nsection); % process data in sections
        iStart           = 1;
        
        
        for iSection = 1:nsection
            data_section    = data(:,iStart:iStart + nsamples_section - 1);
            data_section    = {data_section'};
            spikes          = ss_detect(data_section,spikes);
            iStart          = iStart + nsamples_section;
        end
        
        %clear data data_section
        
        spikes = ss_align(spikes);
        spikes = ss_kmeans(spikes);
        spikes = ss_energy(spikes);
        spikes = ss_aggregate(spikes);
        
        metadata.tetrode = iTetrode;
        save(['Spikes_'  metadata.subject '_' metadata.session '_T' num2str(iTetrode) '_AG' num2str(spikes.params.agg_cutoff) '_ST' num2str(spikes.params.thresh) '.mat'],'spikes','metadata');
    end
    
    toc
    
end

end