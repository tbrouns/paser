function spikes = psr_get_stimulus_window(spikes,params)

tPre  = params.t_win(1)   / 1000;
tPost = params.t_win(end) / 1000;

% Convert
Tmax  = spikes.info.dur;
Fs    = spikes.Fs;
Nmax  = floor(Fs * Tmax) + 1;
sPre  = round(Fs * tPre);  % pre- stimulus window
sPost = round(Fs * tPost); % post-stimulus window

% Extract stimulus windows

spikes.spiketimes = spikes.spiketimes - spikes.info.trialonset;
spikePoints = round(Fs * spikes.spiketimes) + 1;
stimulusTimes = spikes.info.stimtimes{1}(:,1) + spikes.info.stimoffset;
nTrials = length(stimulusTimes); % number of stimulus onsets
    
spikes.trials = zeros(size(spikes.spiketimes));
spikes.info.trialonset = stimulusTimes + tPre;

spikeSignal = zeros(Nmax,1);
spikeSignal(spikePoints) = 1;

if (~isempty(stimulusTimes))
    
    stimulusPoints = round(Fs * stimulusTimes) + 1;
    if (size(stimulusPoints,1) > size(stimulusPoints,2)); stimulusPoints = stimulusPoints'; end
    
    % Find which spikes to remove
    
    stimTrialIDs = bsxfun(@plus,stimulusPoints,(sPre + 1 : sPost)');
    stimTrialIDs(stimTrialIDs < 1)    = 1;
    stimTrialIDs(stimTrialIDs > Nmax) = Nmax;
        
    for iTrial = 1:nTrials
        I = stimTrialIDs(:,iTrial);       
        spikeIDs = extractSpikeIDs(I,spikeSignal);
        spikes.trials(spikeIDs) = iTrial;
    end
    
    del = find(spikes.trials == 0);
    spikes = psr_sst_remove_spikes(spikes,del,'delete');
    
end 

end

function spikeIDs = extractSpikeIDs(I,spikeSignal)

I = I(spikeSignal(I) == 1);
spikeSignal(I) = spikeSignal(I) + 1;
spikeIDs = spikeSignal(spikeSignal > 0) == 2;

end