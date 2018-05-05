function psr_ms_denoise_spk(filesSaved,filesTemp,parameters)

% Spike removal parameters

spikes  = [];
thresh  = parameters.ms.denoise.spk.thresh;
nProbes = size(filesTemp,1);
nTrials = size(filesTemp,2);

% Load info from first file

filename = filesSaved{1};
load(filename,'spikes');

Tmax  = sum(spikes.info.dur);
Fs    = spikes.Fs;
Nmax  = floor(Fs * Tmax) + 1;
sBin  = round(Fs * parameters.ms.denoise.spk.tbin);
sWin  = round(Fs * parameters.ms.denoise.spk.twin);
sBinHalf = round(0.5 * sBin);
trialonsets = [0;cumsum(spikes.info.dur)];

% Convert stimulus timings

I = [];
stimulusPointsAll = cell(0);
for iTrial = 1:nTrials
    load(filesTemp{1,iTrial},'metadata');
    if (metadata.stimulus > 0)
        stimulusPointsAll{end+1} = metadata.stimtimes{1} + trialonsets(iTrial);
        I(end+1) = metadata.stimulus >= parameters.ms.denoise.spk.stim;
    end
end
I = logical(I);

stimulusPoints    = stimulusPointsAll(I);
stimulusPointsAll = cat(1, stimulusPointsAll{:});
stimulusPoints    = cat(1, stimulusPoints{:});
stimulusPointsAll = stimulusPointsAll(:,1);
stimulusPoints    = stimulusPoints(:,1);

% Convert to data point index

stimulusPointsAll = round(Fs * stimulusPointsAll) + 1;
stimulusPoints    = round(Fs * stimulusPoints)    + 1;

stimTrialIDs = bsxfun(@plus,stimulusPoints,-sWin:sWin);
stimTrialIDs(stimTrialIDs < 1)    = 1;
stimTrialIDs(stimTrialIDs > Nmax) = Nmax;

% Load all spikes in vector
spikeSignal = zeros(Nmax,1); % Binary vector of spike points
for iProbe = 1:nProbes
    filename = filesSaved{iProbe};
    load(filename);
    spikePoints = round(Fs * spikes.spiketimes) + 1;
    spikePoints(spikePoints > Nmax) = [];
    spikeSignal(spikePoints) = spikeSignal(spikePoints) + 1;
end

% Set parameters

nSpikes = sum(spikeSignal);
nExpect = sBin * nSpikes / Nmax; % Expected number of spikes in a bin
thresh  = thresh * nExpect;

% Find which spikes to remove

spikeSignalStim = spikeSignal(stimTrialIDs);
spikeSignalStim = sum(spikeSignalStim,1) / size(stimTrialIDs,1);
u = ones(1,sBin);
spikeSignalConv = conv(full(spikeSignalStim)',u,'same');
indicesThresh = find(spikeSignalConv > thresh); % Indices of above-threshold firing rates in window
if (~isempty(indicesThresh))
    
    indicesThresh = bsxfun(@plus,indicesThresh,-sBinHalf:sBinHalf);
    indicesThresh = indicesThresh(:);
    indicesThresh = unique(indicesThresh);
    
    indicesDelete = bsxfun(@plus,stimulusPointsAll,(indicesThresh - sWin)');
    indicesDelete = indicesDelete(:);
    
    % Delete spikes
        
    for iProbe = 1:nProbes
        
        filename = filesSaved{iProbe};
        load(filename);
        
        spikePoints = round(Fs * spikes.spiketimes) + 1;
        spikeSignal = zeros(Nmax,1); % Binary vector of spike points
        spikeSignal(spikePoints) = 1;
        
        del = indicesDelete(spikeSignal(indicesDelete) == 1);
        spikeSignal(del) = spikeSignal(del) + 1;
        del = spikeSignal(spikeSignal > 0) == 2;
        spikes.delete.mse_cluster = del;
        spikes.info.msa = (indicesThresh - sWin) / Fs;

        save(filename,'spikes','-append'); % Save spikes
    end
end

end