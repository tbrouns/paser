function clusters = psr_sst_cluster_features(spikes,freq,metadata,parameters)

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

clusterMin = parameters.cluster.min_spikes;
threshold  = parameters.cluster.thresh * (spikes.info.thresh / parameters.spikes.thresh);
Fs         = spikes.Fs;

% Extract local field potential artifacts

artifacts = [];
if (isfield(freq,'artifacts'))
    ntrials = size(freq,2);
    for iTrial = 1:ntrials
        artifacts = [artifacts;freq(iTrial).artifacts + metadata.trialonset(iTrial)]; %#ok
    end
    artifacts = getArtifactVector(spikes,artifacts);
end

% Convert to single

if (isa(spikes.waveforms,'int16'))
    spikes.waveforms = psr_single(spikes.waveforms,parameters);
end

% Filter spikes

if (parameters.cluster.remove_corr); spikes = psr_sst_filter_corr(spikes,parameters,'array'); end
if (parameters.cluster.remove_amp);  spikes = psr_sst_filter_amp (spikes,parameters,'array'); end
if (parameters.cluster.remove_rpv);  spikes = psr_sst_filter_rpv (spikes,parameters,'array');
else,                                spikes = psr_sst_filter_rpv (spikes,parameters, 'none');
end

if (isfield(spikes,'removed'));      spikes = psr_sst_spike_removal(spikes,find(spikes.removed),'delete'); end

% Calculate features of clusters

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);

clusters = [];
clusters.bhattacharyya = NaN(nClust,nClust);
clusters.mahal         = NaN(nClust,nClust);

for iClust = 1:nClust
    
    % Extract cluster ID
    
    clusterID = clusterIDs(iClust);
    nspikes   = sum(spikes.assigns == clusterID);
    
    % Calculate individual cluster metrics
    
    [p,~,~,~,~]  = ss_undetected        (spikes, clusterID, parameters);             
    co           = ss_censored          (spikes, clusterID, parameters); 
    [xc,lagMax]  = crossCorrelation     (spikes, clusterID, parameters); 
    [k,M,m]      = checkThreshold       (spikes, clusterID, threshold);
    f            = getFiringRate        (spikes, clusterID);
    od           = getArtifactOverlap   (spikes, clusterID, artifacts); 
    snr          = getSignalToNoise     (spikes, clusterID);
    [auc,N]  	 = fitPoisson           (spikes, clusterID, parameters); 
    corrGlobal   = getCorrelationGlobal (spikes, clusterID); 
    corrProbe    = getCorrelationProbe  (spikes, clusterID);
    
    clusters.vars(iClust).id       = clusterID;
    clusters.vars(iClust).rpv      = spikes.rpvs(iClust);
    clusters.vars(iClust).missing  = p;
    clusters.vars(iClust).co       = co;
    clusters.vars(iClust).nspikes  = nspikes;
    clusters.vars(iClust).fspikes  = nspikes / length(spikes.spiketimes);
    clusters.vars(iClust).xc       = xc;           
    clusters.vars(iClust).xc_lag   = 1000 * (lagMax / Fs);
    clusters.vars(iClust).p_auc    = auc;
    clusters.vars(iClust).p_dist   = N;
    clusters.vars(iClust).amp      = M;
    clusters.vars(iClust).amp_rel  = m;
    clusters.vars(iClust).chans    = k;
    clusters.vars(iClust).frate    = f;
    clusters.vars(iClust).artifact = od;
    clusters.vars(iClust).corrG    = corrGlobal;
    clusters.vars(iClust).corrP    = corrProbe;
    clusters.vars(iClust).snr      = snr;
    
    % Pair-wise cluster metrics (unreliable)
    
    PC_1 = psr_pca(spikes,parameters.cluster.pca_dims,find(spikes.assigns == clusterID));
    for jcluster = (iClust+1):nClust
        PC_2 = psr_pca(spikes,parameters.cluster.pca_dims,find(spikes.assigns == clusterIDs(jcluster)));
        if (length(PC_1) > parameters.cluster.pca_dims * clusterMin && length(PC_2) > parameters.cluster.pca_dims * clusterMin)
            M = mahal(PC_1,PC_2);
            B = bhattacharyya(PC_1,PC_2);
            clusters.mahal(iClust,jcluster)         = mean(M(:));
            clusters.bhattacharyya(iClust,jcluster) = B;            
        end
    end
end

% Mixture of drifting t-distribution model for sorting spikes and measuring unit isolation

metrics  = psr_sst_cluster_MDT(spikes,parameters);
ID1      = [metrics.id];
ID2      = [clusters.vars.id];
nClusts  = length(ID1);

for iClust = 1:nClusts
    I = find(ID2 == ID1(iClust));
    if (~isempty(I))
        clusters.vars(I).Lratio = metrics(iClust).Lratio;
        clusters.vars(I).IsoDis = metrics(iClust).IsoDis;
        clusters.vars(I).FP_t   = metrics(iClust).FP_t;
        clusters.vars(I).FN_t   = metrics(iClust).FN_t;
        clusters.vars(I).FP_g   = metrics(iClust).FP_g;
        clusters.vars(I).FN_g   = metrics(iClust).FN_g;
    end
end

end

function RPVtotal = getRPVtotal(spikes,show)

which = ismember(spikes.assigns, show);
RPVtotal = sum(spikes.rpvs(which));

end

function [xc,lagMax] = crossCorrelation(spikes, show, parameters)

clus      = ismember(spikes.assigns, show);
nchan     = size(spikes.waveforms,3);
nsamples  = size(spikes.waveforms,2);
threshold = parameters.cluster.thresh_xcorr * (spikes.info.thresh / parameters.spikes.thresh);

%% Calculate cross-correlation between consecutive spike channels

waves = spikes.waveforms(clus,:,:);
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

function [auc,N] = fitPoisson(spikes, show, parameters)

show = ismember(spikes.assigns, show);

% Parse inputs

spiketimes = sort(spikes.spiketimes(show));
nspikes = length(spiketimes);
twin    = parameters.cluster.stability_win;
tlength = spikes.info.dur;
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

kmax = max(counts);
edges = -0.5:kmax+0.5;
N = histcounts(counts, edges, 'Normalization', 'probability');

% Poisson distribution

lambda = frate * twin;
xmax = max(counts);
x = 0:xmax;
y = poisspdf(x,lambda);

% Find difference between theoretical and empirical Poisson distributions

A = N;
A(A > y) = y(A > y); 
auc = sum(A) / sum(N); % area under Poisson curve (relative to total area)

end

function [k,M,m] = checkThreshold(spikes, show, threshold)

signThresh   = sign(mean(threshold));
clus         = ismember(spikes.assigns, show);
memberwaves  = spikes.waveforms(clus,:,:);
memberwaves  = signThresh * squeeze(mean(memberwaves,1));
M            = max(memberwaves); % maximum amplitude per channel
k            = find(M > abs(threshold)); % channels that cross threshold with mean amplitude
m            = max(M ./ abs(threshold));
M            = max(M); % maximum mean channel amplitude

end

function f = getFiringRate(spikes,show)

clus       = ismember(spikes.assigns, show);
spiketimes = spikes.spiketimes(clus);
nspikes    = length(spiketimes);
dur        = spikes.info.dur;
f          = nspikes / dur; % in Hz

end

function artifactVector = getArtifactVector(spikes,artifacts)

Fs        = spikes.Fs;
duration  = spikes.info.dur;
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

function overlap_diff = getArtifactOverlap(spikes,show,artifacts)

overlap_diff = 0;
if (~isempty(artifacts))
    
    % Finds the degree to which spikes overlap with artifacts in LFP signal
    
    Fs         = spikes.Fs;
    which      = ismember(spikes.assigns, show);
    spiketimes = spikes.spiketimes(which);
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

end

function snr = getSignalToNoise(spikes,show)

% See: Reliability of Signals From a Chronically Implanted Silicon-Based
% Electrode Array in Non-Human Primate Primary Motor Cortex (2005)

which    = ismember(spikes.assigns, show);
waves    = spikes.waveforms(which,:);
wavesAvg = mean(waves);
A        = max(wavesAvg) - min(wavesAvg);
noise    = bsxfun(@minus,waves,wavesAvg);
noise    = noise(:);
noise    = 2 * std(noise);
snr      = A / noise;

end

function corrGlobal = getCorrelationGlobal(spikes,show)
    which = ismember(spikes.assigns, show);
    corrGlobal = nanmean(spikes.correlations(which));
end

function corrProbe = getCorrelationProbe(spikes,show)
    nChans = size(spikes.waveforms,3);
    which = ismember(spikes.assigns, show);
    waves = spikes.waveforms(which,:,:);
    waves = permute(waves,[2 1 3]);
    waves = reshape(waves,[],nChans);    
    R = triu(corr(waves),1);
    R = R(triu(true(size(R)),1));
    corrProbe = min(R);
end

%------------- END OF CODE --------------