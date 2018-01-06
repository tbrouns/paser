function metrics = psr_sst_cluster_features(spikes,freq,parameters)

% PSR_SST_CLUSTER_FEATURES - Calculates quality metrics for spike clusters.
% This function calcules a number of statistics to assess cluster quality
% for every spike cluster obtained with spike sorting.
%
% Syntax: clusters = psr_sst_cluster_features(spikes,freq,metadata,parameters)
%
% Inputs:
%    spikes     - See README
%    freq       - See README
%    metadata   - See README
%    parameters - See README
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
threshold  = parameters.cluster.thresh * (mean(spikes.info.thresh) / parameters.spikes.thresh);
Fs         = spikes.Fs;
clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);

% Extract local field potential artifacts
artifacts = [];
if (isfield(freq,'artifacts'))
    onsets  = [0;cumsum(spikes.info.dur)];
    ntrials = size(freq,2);
    for iTrial = 1:ntrials
        artifacts = [artifacts;freq(iTrial).artifacts + onsets(iTrial)]; %#ok
    end
    artifacts = getArtifactVector(spikes,artifacts);
end

% Filter spikes
if (isfield(spikes,'delete')); spikes = psr_sst_filter_spikes(spikes,parameters,'delete'); end

% Convert to single
if (isa(spikes.waveforms,'int16'))
    spikes.waveforms = psr_single(spikes.waveforms,parameters); 
end

% Calculate features of clusters

metrics = [];

for iClust = nClust:-1:1
    
    % Extract cluster ID
    
    clusterID = clusterIDs(iClust);
    nspikes   = sum(spikes.assigns == clusterID);
    
    metrics(iClust).id       = clusterID;
    metrics(iClust).nspikes  = nspikes;
    metrics(iClust).fspikes  = nspikes / length(spikes.spiketimes);
    
    if (nspikes < parameters.cluster.min_spikes); continue; end
    
    % Calculate individual cluster metrics
    
    fRPV                  = getRPVs             (spikes, clusterID, parameters);
    [fSub,~,~]            = psr_sst_amp_gaussfit(spikes, clusterID, parameters);
    fCoinciding           = ss_censored         (spikes, clusterID, parameters);
    [xcross,lagMax]       = crossCorrelation    (spikes, clusterID, parameters);
    [pAuc,pDist]          = fitPoisson          (spikes, clusterID, parameters);
    [k,ampAbs,ampRel,p2p] = checkThreshold      (spikes, clusterID, threshold);
    overlapFactor         = getArtifactOverlap  (spikes, clusterID, artifacts);
    firingRate            = getFiringRate       (spikes, clusterID);
    signal2noise          = getSignalNoiseRatio (spikes, clusterID);
    corrGlobal            = getCorrelationGlobal(spikes, clusterID);
    corrProbe             = getCorrelationProbe (spikes, clusterID);
    
    metrics(iClust).rpv      = fRPV;
    metrics(iClust).sub      = fSub;
    metrics(iClust).co       = fCoinciding;
    metrics(iClust).xc       = xcross;
    metrics(iClust).xc_lag   = 1000 * (lagMax / Fs); % [ms]
    metrics(iClust).p_auc    = pAuc;
    metrics(iClust).p_dist   = pDist;
    metrics(iClust).amp      = ampAbs;
    metrics(iClust).amp_rel  = ampRel;
    metrics(iClust).p2p      = p2p;
    metrics(iClust).chans    = k;
    metrics(iClust).frate    = firingRate;
    metrics(iClust).artifact = overlapFactor;
    metrics(iClust).corrG    = corrGlobal;
    metrics(iClust).corrP    = corrProbe;
    metrics(iClust).snr      = signal2noise;
end

% Mixture of drifting t-distribution model for sorting spikes and measuring unit isolation

metricsMDT = psr_sst_cluster_MDT(spikes,parameters);
ID1      = [metricsMDT.id];
ID2      = [metrics.id];
nClusts  = length(ID1);

for iClust = nClusts:-1:1
    I = find(ID2 == ID1(iClust));
    if (~isempty(I))
        metrics(I).Lratio = metricsMDT(iClust).Lratio;
        metrics(I).IsoDis = metricsMDT(iClust).IsoDis;
        metrics(I).FP_t   = metricsMDT(iClust).FP_t;
        metrics(I).FN_t   = metricsMDT(iClust).FN_t;
        metrics(I).FP_g   = metricsMDT(iClust).FP_g;
        metrics(I).FN_g   = metricsMDT(iClust).FN_g;
    end
end

end

function fRPV = getRPVs(spikes, clustID, parameters)

if (isfield(spikes.clusters,'rpvs'))
    clustIDs = spikes.clusters.rpvs(:,2);
    I = find(clustIDs == clustID,1);
    fRPV = spikes.clusters.rpvs(I,1);
else
    spikeIDs   = ismember(spikes.assigns, clustID);
    spiketimes = spikes.spiketimes(spikeIDs);
    rpvs       = diff(spiketimes) <= 0.001 * parameters.spikes.ref_period;
    fRPV       = sum(rpvs) / length(spiketimes);
end


end

function [xc,lagMax] = crossCorrelation(spikes, clustID, parameters)

spikeIDs  = ismember(spikes.assigns, clustID);
nchan     = size(spikes.waveforms,3);
nsamples  = size(spikes.waveforms,2);
threshold = parameters.cluster.thresh_xcorr * (mean(spikes.info.thresh) / parameters.spikes.thresh);

%% Calculate cross-correlation between consecutive spike channels

waves = spikes.waveforms(spikeIDs,:,:);
waves_mean = squeeze(mean(waves,1));
waves = reshape(permute(waves,[2 1 3]),size(waves,1)*size(waves,2),size(waves,3));

chanIDs = zeros(1,nchan);
waves_max = max(abs(waves_mean));
threshold = abs(threshold); % consider both positive and negative peaks
for ichan = 1:nchan
    chanIDs(ichan) = waves_max(ichan) >= threshold(ichan);
end

maxlag = round(0.75 * nsamples);

N  = (nchan / 2) * (nchan - 1);
XC = zeros(2 * maxlag + 1, N);
accept = zeros(1,N);

n  = 1;
for ichan = 1:nchan
    x = waves(:,ichan);
    for jchan = ichan+1:nchan
        y = waves(:,jchan);
        [XC(:,n),lags] = xcorr(x,y,maxlag,'coeff');
        accept(n) = chanIDs(ichan) & chanIDs(jchan);
        n = n + 1;
    end
end

% Ignore low amplitude channels
XC = XC(:,logical(accept));
ncombi = size(XC,2);

% Find location of xcorr peak with largest lag

lagAll = zeros(ncombi,1);

for icombi = 1:ncombi
    xc = XC(:,icombi);
    [pks,loc] = findpeaks(xc);
    if (isempty(pks)); [~,I] = max(xc);
    else,              [~,I] = max(pks); I = loc(I);
    end
    lagAll(icombi) = lags(I);
end

[lagMax,I] = max(lagAll);
if (isempty(lagMax)); lagMax = 0; end
xc = XC(:,I);

end

function [pAuc,pDist] = fitPoisson(spikes, clustID, parameters)

spikeIDs = ismember(spikes.assigns, clustID);

% Parse inputs

spiketimes = sort(spikes.spiketimes(spikeIDs));
nspikes = length(spiketimes);
twin    = parameters.cluster.stability_win;
tlength = sum(spikes.info.dur);
frate   = nspikes / tlength;
nWin    = floor(2 * tlength / twin); % number of windows
counts  = zeros(nWin,1);
tstart  = 1;

% Acquire spike counts
for i = 1:nWin
    tend = tstart + twin; % Update frequency window
    if (tend > tlength) % Check if end of signal is reached
        tstart = tlength - twin;
        tend   = tlength;
    end
    counts(i) = sum(spiketimes >= tstart & spiketimes < tend);
    tstart = tstart + round(0.5 * twin); % Update starting frequency
end

kmax  = max(counts);
edges = -0.5:kmax+0.5;
pDist = histcounts(counts, edges, 'Normalization', 'probability');

% Poisson distribution

lambda = frate * twin;
xmax = max(counts);
x = 0:xmax;
y = poisspdf(x,lambda);

% Find difference between theoretical and empirical Poisson distributions

A = pDist;
A(A > y) = y(A > y);
pAuc = sum(A) / sum(pDist); % area under Poisson curve (relative to total area)

end

function [chanIDs,maxAmp,relAmp,p2p] = checkThreshold(spikes, clustID, threshold)

signThresh = sign(mean(threshold));
spikeIDs   = ismember(spikes.assigns, clustID);
waves      = spikes.waveforms(spikeIDs,:,:);
waves      = signThresh * squeeze(mean(waves,1));
maxAmp     = max(waves); % maximum amplitude per channel
chanIDs    = find(maxAmp > abs(threshold)); % channels that cross threshold with mean amplitude
relAmp     = max(maxAmp ./ abs(threshold));
maxAmp     = max(maxAmp); % maximum mean channel amplitude
p2p        = maxAmp - min(waves(:));

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
nlength    = length(artifacts);

spiketimes(spiketimes < 1)       = [];
spiketimes(spiketimes > nlength) = [];

nspikes = length(spiketimes);

n = sum(artifacts);
N = sum(artifacts(spiketimes));
overlap_expected = n / nlength; % expected overlap in case of uniform spike distribution
overlap_exact    = N / nspikes; % actual fraction of spikes that overlap with artifact

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

function corrGlobal = getCorrelationGlobal(spikes,clustID)

spikeIDs = ismember(spikes.assigns, clustID);
corrGlobal = nanmean(spikes.corr_global(spikeIDs));

end

function corrProbe = getCorrelationProbe(spikes,clustID)

nChans = size(spikes.waveforms,3);
spikeIDs = ismember(spikes.assigns, clustID);
waves = spikes.waveforms(spikeIDs,:,:);
waves = permute(waves,[2 1 3]);
waves = reshape(waves,[],nChans);
R = triu(corr(waves),1);
R = R(triu(true(size(R)),1));
corrProbe = min(R);

end

%------------- END OF CODE --------------