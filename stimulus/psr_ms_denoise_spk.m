function psr_ms_denoise_spk(filesSaved)

% Spike removal parameters

tWin       = 0.050; % Window size [ms]
tBin       = 0.001; % Bin size [ms]
threshStim = 0;    % [V]
threshMin  = 0.25;  % At least X of the trials should have spike in window
thresh     = 8;     % Spike count threshold within window
corrThresh = 0.5;

nProbes = length(filesSaved);
spikes = [];

% Load info from first file

filename = filesSaved{1};
load(filename);

Tmax  = sum(spikes.info.dur);
Fs    = spikes.Fs;
Nmax  = floor(Fs * Tmax) + 1;
sBin  = round(Fs * tBin);
sWin  = round(Fs * tWin);

sBinHalf = round(0.5 * sBin);

stimuli  = cell2mat(metadata.stimulus);
trialIDs_0 = stimuli > 0;

stimuli = stimuli(trialIDs_0);
trialIDs_1 = stimuli > threshStim;

% Convert stimulus timings

itr = 1;
stimulusPointsAll = metadata.stimtimes(trialIDs_0,1);
for iTrial = find(trialIDs_0)'
    stimulusPointsAll{itr} = stimulusPointsAll{itr} + metadata.trialonset(iTrial);
    itr = itr + 1;
end

stimulusPoints    = stimulusPointsAll(trialIDs_1);
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

% Binary vector of spike points

spikeSignal = zeros(Nmax,1);

% Load all spikes in vector
for iProbe = 1:nProbes
    filename = filesSaved{iProbe};
    load(filename);
    spikePoints = round(Fs * spikes.spiketimes) + 1;
    spikeSignal(spikePoints) = spikeSignal(spikePoints) + 1;
end

% Set parameters

nSpikes = length(sum(spikeSignal));
nExpect = sBin * nSpikes / Nmax; % Expected number of spikes in a bin
thresh  = thresh * nExpect;
if (thresh < threshMin); thresh = threshMin; end

% Find which spikes to remove

spikeSignalStim = spikeSignal(stimTrialIDs);
spikeSignalStim = sum(spikeSignalStim,1) / size(stimTrialIDs,1);
u = ones(1,sBin);
spikeSignalConv = conv(full(spikeSignalStim)',u,'same');
indicesThresh = find(spikeSignalConv > thresh); % Indices of above-threshold firing rates in window
if (~isempty(indicesThresh))
    dIndicesThresh = diff(indicesThresh);
    indicesSplit  = [0,find(dIndicesThresh > 1),length(dIndicesThresh)];
    indicesCentre = [];
    for i = 1:length(indicesSplit)-1
        indicesAdjacent = indicesSplit(i)+1:indicesSplit(i+1);
        ConvAdjacent = spikeSignalConv(indicesThresh(indicesAdjacent));
        SUCCESS = false;
        if (length(ConvAdjacent) >= 3)
            [~,locs] = findpeaks(ConvAdjacent);
            if (~isempty(locs)); SUCCESS = true; end
        end
        if (~SUCCESS); locs = round(0.5 * length(ConvAdjacent)); end
        indicesCentre = [indicesCentre;indicesThresh(locs)];
    end
    
    dI = diff(indicesCentre);
    indicesSplit = [0,find(dI > sBinHalf),length(dI)];
    
    % Remove spikes for each probe
    for iProbe = 1:nProbes
        
        filename = filesSaved{iProbe};
        load(filename);
        spikePoints = round(Fs * spikes.spiketimes) + 1;
        
        for iSplit = 1:length(indicesSplit)-1
            itr = 1;
            I  = indicesSplit(iSplit)+1:indicesSplit(iSplit+1);
            n1 = length(indicesCentre(I));
            n2 = length(stimulusPointsAll);
            indicesDelete = zeros(n1*n2,1);
            for i = 1:length(indicesCentre)
                I = itr : itr + n2 - 1;
                indicesDelete(I,1) = stimulusPointsAll - sWin + indicesCentre(i);
                itr = itr + n2;
            end
            
            IDs = bsxfun(@plus,indicesDelete,-sBinHalf:sBinHalf);
            
            spikeSignal = zeros(Nmax,1);
            spikeSignal(spikePoints) = 1;
            IDs = IDs(spikeSignal(IDs) == 1);
            spikeSignal(IDs) = spikeSignal(IDs) + 1;
            IDs = spikeSignal(spikeSignal > 0) == 2;
            spikeIDs = find(IDs);
            
            waveforms = psr_int16_to_single(spikes.waveforms(spikeIDs,:),parameters);
            waveforms = waveforms';
            waveformMedian = median(waveforms,2);
            
            nspikes = size(waveforms,1);
            R = zeros(nspikes,1);
            for iSpike = 1:nspikes
                R(iSpike) = corr(waveformMedian,waveforms);
            end
            
            spikeIDsAll = [spikeIDsAll;spikeIDs(R > corrThresh)];

        end
                    
        spikes = psr_sst_remove_spikes(spikes,spikeIDsAll,'delete');
        spikes.info.msa = (indicesCentre - sWin) / Fs;
                
        save(filename,'spikes','-append'); % Save spikes
    end
end

end