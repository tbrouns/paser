function dataFT = psr_convert2fieldtrip(data,timestamps,stimTimes,parameters)

Fr      = parameters.Fr;
nchans  = size(data,1);
labels  = 1:nchans;

method    = stimTimes{2};
stimTimes = stimTimes{1};

% Set onset and offset for every trial

if (strcmp(method,'onset')) % only use onset
    
    onset  = parameters.lfp.trial.onset  - parameters.lfp.trial.padding;
    offset = parameters.lfp.trial.offset + parameters.lfp.trial.padding;
        
    onset     = round(Fr * onset);
    offset    = round(Fr * offset);
    stimTimes = round(Fr * stimTimes(:,1));
    
    stimOnset  = stimTimes + onset;
    stimOffset = stimTimes + offset;    
        
elseif (strcmp(method,'interval')) % Use both onset and offset to specify interval
        
    onset = 0;
    stimOnset  = round(Fr * stimTimes(:,1));
    stimOffset = round(Fr * stimTimes(:,2));
    
end

% Ignore trials outside boundary
nlength = round(Fr * timestamps(end));
del = stimOnset < 1 | stimOffset > nlength;
stimOnset (del) = [];
stimOffset(del) = [];
    
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
    cfg.offset = onset;
    dataFT     = ft_redefinetrial(cfg,dataFT);
    
    % Set NaNs to zero
    
    nTrials = size(dataFT.trial,2);
    for iTrial = 1:nTrials
        data_trial = sum(dataFT.trial{iTrial});
        dataFT.trial{iTrial}(:,isnan(data_trial)) = 0;
    end
    
end