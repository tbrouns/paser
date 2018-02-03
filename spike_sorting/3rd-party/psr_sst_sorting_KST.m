function spikes = psr_sst_sorting_KST(data,parameters,savePath)

disp('Running Kilosort...');

addpath(genpath(parameters.path.kst));

% Dummy parameters

nchans      = size(data,1);
chanMap     = 1:nchans;
chanMap0ind = chanMap - 1;    %#ok
connected   = true(nchans,1); %#ok
xcoords     = ones(nchans,1); %#ok
ycoords     = (1:nchans)';    %#ok 
kcoords     = ones(nchans,1); %#ok
fs          = parameters.Fs;  %#ok

save(fullfile(savePath, 'chanMap.mat'), ...
    'chanMap',...
    'chanMap0ind',...
    'connected',...
    'xcoords',...
    'ycoords',...
    'kcoords',...
    'fs');

% Kilosort parameters

ops = parameters.sorting.kst;
ops.Nchan = nchans;
ops = kilosort_config(ops,parameters,savePath);

% Calculate threshold for spike initialization

sd = psr_mad(single(data'));
sd = mean(sd); % Take mean across channels, because we need scalar
ops.scaleproc = sd;

% Convert data to DAT filetype with int16 data

ops.data = data; 
fidout = fopen(ops.fbinary, 'w');
fwrite(fidout, ops.data, 'int16');
fclose(fidout);
ops = rmfield(ops,'data');

% NOTE: we do not use Kilosort's "convertOpenEphysToRawBInary" to load from
% CONTINUOUS files here, because it is incompatible with our data
% processing pipeline.

% Kilosort pipeline

try
    [rez, DATA, uproj] = preprocessData(ops);            % preprocess data and extract spikes for initialization
    rez                = fitTemplates(rez, DATA, uproj); % fit templates iteratively
    rez                = fullMPMU(rez, DATA);            % extract final spike times (overlapping extraction)
catch
    disp('Kilosort error. No spikes detected.');
    rez = [];
end

if (~isempty(rez))
    % extract info from rez
    spiketimes = rez.st3(:,1);
    if (size(rez.st3,2) >= 5); assigns = 1 + rez.st3(:,5);
    else,                      assigns =     rez.st3(:,2);
    end
    spikes = psr_convert2spikes(data,spiketimes,assigns,parameters);
    rez = rmfield(rez,'st3'); % Remove now-redundant field
    spikes.info.kst = rez;
else               
    spikes = [];
end

% Clean-up
fclose('all'); 
delete(ops.fbinary);
delete(ops.fproc);
delete(ops.chanMap);
rmpath(genpath(parameters.path.kst));

disp('Kilosort completed');

end

function ops = kilosort_config(ops,parameters,savePath)

ops.parfor      = 0; % whether to use parfor to accelerate some parts of the algorithm
ops.showfigures = 0; % whether to plot figures during optimization

ops.datatype = 'dat';  % binary ('dat', 'bin') or 'openEphys'
ops.fbinary  = fullfile(savePath, 'raw_data.dat'); % will be created for 'openEphys'
ops.fproc    = fullfile(savePath, 'temp_wh.dat'); % residual from RAM of preprocessed data

% define the channel map as a filename (string) or simply an array
ops.chanMap  = fullfile(savePath, 'chanMap.mat'); % make this file using createChannelMapFile.m

ops.nNeighPC = ops.Nchan; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)
ops.nNeigh   = 16;        % visualization only (Phy): number of neighboring templates to retain projections of (16)

ops.fshigh = parameters.spikes.bp.lower; % frequency for high pass filtering
ops.fslow  = parameters.spikes.bp.upper; % frequency for low pass filtering (optional)

% options for channel whitening
ops.whitening      = 'noSpikes'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
ops.nSkipCov       = 1;          % compute whitening matrix from every N-th batch (1)
ops.whiteningRange = Inf;        % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)

ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info).

% load predefined principal components (visualization only (Phy): used for features)
dd       = load('PCspikes2.mat'); % you might want to recompute this from your own data
ops.wPCA = dd.Wi(:,1:7);   % PCs

% options for posthoc merges (under construction)
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
ops.epu     = Inf;

ops.ForceMaxRAMforDat = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.

end