function dataFT = psr_convert2fieldtrip(data,timestamps,stimTimes,parameters)

Fr      = parameters.Fr;
nchans  = size(data,1);
labels  = 1:nchans;

% Set onset and offset for every trial

if (size(stimTimes,2) == 1) % Only onset specified
    onset  = parameters.lfp.trial_onset  - parameters.lfp.trial_padding;
    offset = parameters.lfp.trial_offset + parameters.lfp.trial_padding;
    stimOnset  = stimTimes + onset;
    stimOffset = stimTimes + offset;
else % Both onset and offset specified
    onset = 0;
    stimOnset  = stimTimes(:,1);
    stimOffset = stimTimes(:,2);
end

offset = round(onset  * Fr); % offset for trial timescale

% Convert to sample number
stimOnset  = round(stimOnset  * Fr);
stimOffset = round(stimOffset * Fr);
nlength = round(Fr * timestamps(end));
stimOnset (stimOnset  < 1) = 1;
stimOffset(stimOffset > nlength) = nlength;
    
% Function that converts dataformat to field trip format
dataFT         = [];
dataFT.label   = strtrim(cellstr(num2str(labels'))'); % cell-array containing strings, Nchan*1
dataFT.fsample = Fr;           % sampling frequency in Hz, single number
dataFT.trial   = {data};       % cell-array containing a data matrix for each 
                               % trial (1 X Ntrial), each data matrix is a Nchan*Nsamples matrix 
dataFT.time    = {timestamps}; % cell-array containing a time axis for each 
                               % trial (1 X Ntrial), each time axis is a 1*Nsamples vector          
dataFT.sampleinfo = [1 size(data,2)]; % array containing [startsample endsample] of data

if (~isempty(stimOnset))
    
    cfg = [];
    cfg.trl = [stimOnset,stimOffset,zeros(size(stimOnset))];
    dataFT  = ft_redefinetrial(cfg,dataFT);
    
    cfg = [];
    cfg.offset = offset;
    dataFT     = ft_redefinetrial(cfg,dataFT);
    
    % Set NaNs to zero
    
    ntrials = size(dataFT.trial,2);
    for iTrial = 1:ntrials
        data_trial = sum(dataFT.trial{iTrial});
        dataFT.trial{iTrial}(:,isnan(data_trial)) = 0;
    end 
    
end