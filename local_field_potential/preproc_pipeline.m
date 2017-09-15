function preproc_pipeline(fpath)

if nargin < 1; fpath = []; end

close all

pattern = 'CH';
ext     = '.continuous';

% Find files

files_unsorted = dir([fpath '*' pattern '*' ext]);
files_unsorted = char(files_unsorted.name);

% sort files

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
 
% filter

% per tetrode

numtets = numfiles / 4;

for iTetrode = 1:numtets;        
    for iElectrode = 1:4;
        
        iFile = (iTetrode - 1) * 4 + iElectrode;
        
        filename = [fpath files{iFile}];
        filename = strtrim(filename);
        
        [data, timestamps, info] = load_open_ephys_data(filename);
        
        Fsample  = info.header.sampleRate;
        
%         nsamples = t_max * Fsample;
%         data        = data      (1:nsamples);
%         timestamps  = timestamps(1:nsamples);
        
        % data in microvolts
        % timestamps in seconds
        
        timestamps  = timestamps - timestamps(1);
        
        if (iElectrode == 1) % new tetrode
            labels       = cell(4,1);
            data_tetrode = zeros(4,length(data));
        end
        
        labels{iElectrode}         = num2str(iFile);
        data_tetrode(iElectrode,:) = data;
    
        data = data';               %#ok
        sr   = Fsample;             %#ok
     
%         save(['data_electrode_' sprintf('%04d',iFile) '.mat'],'data','sr', '-v7.3');
    
    end
        
    data_input            = [];
    data_input.label      = labels;
    data_input.fsample    = Fsample;
    data_input.trial      = {data_tetrode};
    data_input.time       = {timestamps'};
    data_input.sampleinfo = [1 length(data_tetrode(1,:))];
        
    data = data_input.trial{1}; %#ok
    sr   = Fsample;             %#ok
    
%     save(['data_tetrode_' sprintf('%04d',iTetrode) '.mat'], 'data',       'sr'   , '-v7.3');
    save(['raw_FT_'   sprintf('%04d',iTetrode) '.mat'], 'data_input'         , '-v7.3');
            
end

end