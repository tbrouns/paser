function rez = ept_sst_sorting_KST(parameters,files,loadPath,savePath)

% KiloSort

addpath(genpath(parameters.sorting.kst.path));

% Dummy parameters

nchans      = size(files,1);
chanMap     = 1:nchans;
chanMap0ind = chanMap - 1;      %#ok
connected   = true(nchans,1);   %#ok
xcoords     = ones(nchans,1);   %#ok
ycoords     = (1:nchans)';      %#ok 
kcoords     = ones(nchans,1);   %#ok
fs          = parameters.Fs;    %#ok

save(fullfile(savePath, 'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs');

ops = kilosort_config(parameters,files,loadPath,savePath);
[rez, DATA, uproj] = preprocessData(ops);            % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj); % fit templates iteratively
rez                = fullMPMU(rez, DATA);            % extract final spike times (overlapping extraction)

% DATA = permute(DATA,[2 1 3]);
% DATA = DATA(:,:);

rmpath(genpath(parameters.sorting.kst.path));

end