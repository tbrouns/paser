function dataFT = ept_convert2fieldtrip(data,timestamps,stimOnset,parameters,Fs)

N       = 2;
nchans  = size(data,1);
tlength = parameters.lfp.trial_length - parameters.lfp.base_onset; 
trial_length = round(N * tlength * Fs); % read in twice as much
labels = 1:nchans;
 
% Function that converts dataformat to field trip format
dataFT = [];
dataFT.label   = strtrim(cellstr(num2str(labels'))'); % cell-array containing strings, Nchan*1
dataFT.fsample = Fs;           % sampling frequency in Hz, single number
dataFT.trial   = {data};       % cell-array containing a data matrix for each 
                               % trial (1 X Ntrial), each data matrix is a Nchan*Nsamples matrix 
dataFT.time    = {timestamps}; % cell-array containing a time axis for each 
                               % trial (1 X Ntrial), each time axis is a 1*Nsamples vector          
dataFT.sampleinfo = [1 size(data,2)]; % array containing [startsample endsample] of data


if (~isempty(stimOnset))

    stimOnset  = (round(stimOnset * Fs))';
    stimOffset = stimOnset + trial_length;
    
    cfg = [];
    cfg.trl    = [stimOnset,stimOffset,zeros(size(stimOnset))]; % MFA start times until next start time
    dataFT = ft_redefinetrial(cfg,dataFT);

    cfg = [];
    cfg.offset = round(N * parameters.lfp.base_onset * Fs);
    dataFT = ft_redefinetrial(cfg,dataFT);

    % Set NaNs to zero
    
    ntrials = size(dataFT.trial,2);
    for iTrial = 1:ntrials
        data_trial = sum(dataFT.trial{iTrial});
        dataFT.trial{iTrial}(:,isnan(data_trial)) = 0;
    end 
    
end