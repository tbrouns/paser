function spikes = psr_ms_denoise_off(spikes,metadata,parameters)

% Load info from first file

Tmax  = sum(spikes.info.dur);
Fs    = spikes.Fs;
Nmax  = floor(Fs * Tmax) + 1;
sWin  = round(Fs * parameters.ms.denoise.off.twin / 1000);
trialonsets = [0;cumsum(spikes.info.dur)];

% Convert stimulus timings

stimulusPoints = cell(0,0);
nTrials = size(metadata.stimtimes,1);
for iTrial = 1:nTrials
    if (metadata.stimulus{iTrial} > 0)
        stimulusPoints{end+1} = metadata.stimtimes{iTrial,1} + trialonsets(iTrial);
    end
end
stimulusPoints = cat(1, stimulusPoints{:});

if (~isempty(stimulusPoints))
    
    stimulusPoints = stimulusPoints(:,1);
    stimulusPoints = round(Fs * stimulusPoints) + 1;
    
    offsets  = round(Fs * metadata.stimoffsetProbe);
    nOffsets = length(offsets);
    
    % Load all spikes in vector
    spikePoints = round(Fs * spikes.spiketimes) + 1;    
    deleteIDs   = false(size(spikePoints));
    
    for iOffset = 1:nOffsets
        
        % Convert to data point index
        
        I = (sWin(1):sWin(2)) + offsets(iOffset);
        offTrialIDs = bsxfun(@plus,stimulusPoints,I);
        offTrialIDs(offTrialIDs < 1)    = 1;
        offTrialIDs(offTrialIDs > Nmax) = Nmax;
        
        deleteIDs = deleteIDs | psr_get_spike_ids(spikePoints,offTrialIDs(:));
        
    end
        
    spikes.delete.off = deleteIDs;
    
end
end