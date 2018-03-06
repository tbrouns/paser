function spikes = psr_sst_cluster_features(spikes,parameters)

% PSR_SST_CLUSTER_FEATURES - Calculates quality metrics for spike clusters.
% This function calcules a number of statistics to assess cluster quality
% for every spike cluster obtained with spike sorting.
%
% Syntax: clusters = psr_sst_cluster_features(spikes,freq,metadata,parameters)
%
% Inputs:
%    spikes     - See README
%    parameters - See README
%    freq       - See README
%
% Output:
%    clusters - See 'spikes.clusters' in README
%
% See also: PSR_SST_CLUSTER_MERGE

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

% Set parameters

Fs         = spikes.Fs;
clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);
spikesOld  = spikes;
artifacts  = [];

% Extract artifacts
if (~psr_isempty_field(spikes,'spikes.info.artifacts'))
    trialOnsets  = [0;cumsum(spikes.info.dur)];
    nTrials = length(spikes.info.artifacts);
    for iTrial = 1:nTrials
        artifactsTemp  = spikes.info.artifacts{iTrial};
        artifacts = [artifacts;artifactsTemp + trialOnsets(iTrial)]; %#ok
    end
    artifacts = getArtifactVector(spikes,artifacts);
end

% Convert to single
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

% Calculate features of clusters
metrics = [];

for iClust = nClust:-1:1
    
    % Extract cluster ID
    
    clusterID = clusterIDs(iClust);
    nspikes   = sum(spikes.assigns == clusterID);
    
    metrics(iClust).id      = clusterID;
    metrics(iClust).nspikes = nspikes;
    metrics(iClust).fspikes = nspikes / length(spikes.spiketimes);
    
    if (nspikes < parameters.cluster.quality.min_spikes); continue; end
    
    % Calculate individual cluster metrics
    
    rpvRatio              = getRPVs                  (spikes, clusterID, parameters);
    [subThreshRatio,~,~]  = psr_sst_amp_gaussfit     (spikes, clusterID, parameters);
    overlapRatio          = psr_sst_spike_overlap    (spikes, clusterID, parameters);
    [xcross,lagMax]       = crossCorrelation         (spikes, clusterID, parameters);
    [cAuc,xDist,yDist]    = psr_sst_cluster_stability(spikes, clusterID, parameters);
    [k,ampAbs,ampRel,p2p] = psr_sst_cluster_amp      (spikes, clusterID, parameters);
    overlapFactor         = getArtifactOverlap       (spikes, clusterID, artifacts);
    firingRate            = getFiringRate            (spikes, clusterID);
    signal2noise          = getSignalNoiseRatio      (spikes, clusterID);
    
    metrics(iClust).frate    = firingRate;
    metrics(iClust).rpv      = rpvRatio;
    metrics(iClust).sub      = subThreshRatio;
    metrics(iClust).co       = overlapRatio;
    metrics(iClust).xc       = xcross;
    metrics(iClust).xcLag    = 1000 * (lagMax / Fs); % [ms]
    metrics(iClust).cAuc     = cAuc;
    metrics(iClust).cxDist   = xDist;
    metrics(iClust).cyDist   = yDist;
    metrics(iClust).amp      = ampAbs;
    metrics(iClust).ampRel   = ampRel;
    metrics(iClust).p2p      = p2p;
    metrics(iClust).chans    = k;
    metrics(iClust).artifact = overlapFactor;
    metrics(iClust).snr      = signal2noise;
end

% Return structure
spikes = spikesOld;
spikes.clusters.metrics = metrics;

end

function fRPV = getRPVs(spikes, clustID, parameters)

spikeIDs   = ismember(spikes.assigns, clustID);
spiketimes = spikes.spiketimes(spikeIDs);
rpvs       = diff(spiketimes) <= 0.001 * parameters.spikes.ref_period;
fRPV       = sum(rpvs) / length(spiketimes);

end

function [crossCorr,lagMax] = crossCorrelation(spikes, clustID, parameters)

spikeIDs  = ismember(spikes.assigns, clustID);
threshold = parameters.cluster.thresh_xcorr * spikes.info.bgn;

%% Calculate cross-correlation between consecutive spike channels

waveforms = spikes.waveforms(spikeIDs,:,:);

nSpikes  = size(waveforms,1);
nSamples = size(waveforms,2);
nChan    = size(waveforms,3);

% Concatenate each spike on every individual channel
chanIDs = zeros(1,nChan);
waveformsMax = max(mean(waveforms,1),[],2);
threshold = abs(threshold); % consider both positive and negative peaks
for iChan = 1:nChan
    chanIDs(iChan) = waveformsMax(iChan) >= threshold(iChan);
end

maxlag = round(0.75 * nSamples);

N = (nChan^2 - nChan) / 2;
crossCorr = zeros(2 * maxlag + 1, N);
accept = false(1,N);

waveforms = reshape(permute(waveforms,[2 1 3]),nSpikes*nSamples,nChan);

itr = 1;
for iChan = 1:nChan
    x = waveforms(:,iChan);
    for jChan = iChan+1:nChan
        y = waveforms(:,jChan);
        [crossCorr(:,itr),lags] = xcorr(x,y,maxlag,'coeff');
        accept(itr) = chanIDs(iChan) & chanIDs(jChan);
        itr = itr + 1;
    end
end

% Ignore low amplitude channels
crossCorr = crossCorr(:,accept);
nCombi = size(crossCorr,2);

% Find location of xcorr peak with largest lag

lagAll = zeros(nCombi,1);

for iCombi = 1:nCombi
    xc = crossCorr(:,iCombi);
    [pks,loc] = findpeaks(xc);
    if (isempty(pks)); [~,I] = max(xc);
    else,              [~,I] = max(pks); I = loc(I);
    end
    lagAll(iCombi) = lags(I);
end

[lagMax,I] = max(lagAll);
if (isempty(lagMax)); lagMax = 0; end
crossCorr = crossCorr(:,I);

end

function f = getFiringRate(spikes,clustID)

spikeIDs   = ismember(spikes.assigns, clustID);
spiketimes = spikes.spiketimes(spikeIDs);
nspikes    = length(spiketimes);
duration   = sum(spikes.info.dur);
f          = nspikes / duration; % in Hz

end

function artifactVector = getArtifactVector(spikes,artifacts)

Fs        = spikes.Fs;
duration  = sum(spikes.info.dur);
nlength   = round(Fs * duration)  + 1;
artifacts = round(Fs * artifacts) + 1;
artifacts(artifacts < 1)       = 1;
artifacts(artifacts > nlength) = nlength;
artifactVector = false(nlength,1);
nArtifacts = size(artifacts,1);
for iArtifact = 1:nArtifacts
    artifactVector(artifacts(iArtifact,1):artifacts(iArtifact,2)) = true;
end
end

function overlap_diff = getArtifactOverlap(spikes,clustID,artifacts)

overlap_diff = 0;
if (isempty(artifacts)); return; end

% Finds the degree to which spikes overlap with artifacts in LFP signal

Fs         = spikes.Fs;
spikeIDs   = ismember(spikes.assigns, clustID);
spiketimes = spikes.spiketimes(spikeIDs);
spiketimes = round(Fs * spiketimes) + 1; % in sample number
nLength    = length(artifacts);

spiketimes(spiketimes < 1)       = [];
spiketimes(spiketimes > nLength) = [];

nSpikes = length(spiketimes);

n = sum(artifacts);
N = sum(artifacts(spiketimes));
overlap_expected = n / nLength; % expected overlap in case of uniform spike distribution
overlap_exact    = N / nSpikes; % actual fraction of spikes that overlap with artifact

overlap_diff = overlap_exact / overlap_expected;

end

function snr = getSignalNoiseRatio(spikes,clustID)

% See: Reliability of Signals From a Chronically Implanted Silicon-Based
% Electrode Array in Non-Human Primate Primary Motor Cortex (2005)

spikeIDs = ismember(spikes.assigns, clustID);
waves    = spikes.waveforms(spikeIDs,:);
wavesAvg = mean(waves);
A        = max(wavesAvg) - min(wavesAvg);
noise    = bsxfun(@minus,waves,wavesAvg);
noise    = noise(:);
noise    = 2 * std(noise);
snr      = A / noise;

end

%------------- END OF CODE --------------