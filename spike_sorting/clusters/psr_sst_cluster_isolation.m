function spikes = psr_sst_cluster_isolation(spikes,parameters,weights)

% Mixture of drifting t-distribution model for sorting spikes and measuring unit isolation
%
% This function is based on 'demo_isolation_metrics' by Kevin Shan [Copyright (c) 2016]
% See LICENSE in MDT folder (or: https://github.com/kqshan/MoDT)
%
% Major edits by Terence Brouns (2017)

if (nargin < 3); weights = ones(size(spikes.spiketimes)); end

spikesOld = spikes;

% Positive non-zero assigns

n = 0;
while (min(unique(spikes.assigns))) <= 0
    spikes.assigns = spikes.assigns + 1;
    n = n + 1;
end

% Ignore clusters that are too small

nClusts = max(spikes.assigns);
if (isempty(nClusts) || nClusts == 0); return; end

clustIDs = 1:nClusts;
nspikes  = zeros(1,nClusts);
for iClust = 1:nClusts
    nspikes(iClust) = sum(spikes.assigns == clustIDs(iClust));
end

I        = nspikes <= 2 * size(spikes.features,1);
del      = clustIDs( I);
clustIDs = clustIDs(~I); % Keep remaining clusters

spikeIDs = ismember(spikes.assigns,del);
spikes = psr_sst_remove_spikes(spikes,find(spikeIDs),'delete');

% Convert to consecutive cluster assigns

itr = 1;
for iClust = clustIDs
    spikeIDs = ismember(spikes.assigns,iClust);
    spikes.assigns(spikeIDs) = itr;
    itr = itr + 1;
end

nClusts = max(spikes.assigns);
if (isempty(nClusts) || nClusts == 0); return; end

% Data conversion

%  New data structure with fields:
%     spk_Y           [D x N] spikes (N spikes in a D-dimensional feature space)
%     spk_t           [N x 1] spike times (ms)
%     spk_clustId     [N x 1] cluster ID this spike is assigned to
%     weighting       [N x 1] suggested weighting for training on a subset

data             = [];
data.spk_Y       = double(spikes.features);
data.spk_t       = double(spikes.spiketimes' * 1000); % convert to ms
data.spk_clustId = double(spikes.assigns');
data.weighting   = double(weights');

if (isempty(data.spk_Y)); return; end

% Run algorithm

% These are the parameters that we recommend in the paper
nu = parameters.cluster.mdt.nu;                               % t-distribution nu parameter (smaller = heavier tails)
q_perhour = parameters.cluster.mdt.q_perhour;                 % Drift regularization (smaller = more smoothing)
timeframe_minutes = parameters.cluster.mdt.timeframe_minutes; % Time frame duration (mostly a computational thing)

% Construct an MoDT object using these parameters
q_perframe = q_perhour * (timeframe_minutes / 60);
model = MoDT('nu',nu,'Q',q_perframe);

% Attach the data to the model
timeframe_ms = 60e3 * timeframe_minutes;
model.attachData(data.spk_Y, data.spk_t,'frameDur',timeframe_ms);

% Fit the model parameters based on our spike assignments
clustAssigned = data.spk_clustId;
MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID); % Disable warning
model.initFromAssign(clustAssigned);
warning('on',  MSGID);

% Obtain the posterior probability that spike n came from cluster k
posterior = model.getValue('posterior');

% Let's also fit a drifting Gaussian model
gaussModel = model.copy();
gaussModel.setParams('nu',Inf);
gaussModel.initFromAssign(clustAssigned);
[gaussPosterior, gaussMahalSq] = gaussModel.getValue('posterior','mahalDist');

% Report some unit isolation metrics
nClust = model.K;
nDims  = model.D;

% Display the results in sorted order
metrics(nClust).id = [];

for k = 1:nClust
    
    % False positive/negative ratios
    is_assigned_to_k = (clustAssigned == k);
    N_k = sum(is_assigned_to_k);
    otherClusts = [1:k-1, k+1:nClust];
    
    % T-distribution
    prob_came_from_k = posterior(:,k);
    prob_came_from_other = sum(posterior(:,otherClusts), 2);
    falsePosMODT = sum(prob_came_from_other( is_assigned_to_k)) / N_k;
    falseNegMODT = sum(prob_came_from_k    (~is_assigned_to_k)) / N_k;
    
    % Repeat this for the Gaussian
    prob_came_from_k = gaussPosterior(:,k);
    prob_came_from_other = sum(gaussPosterior(:,otherClusts), 2);
    falsePosGauss = sum(prob_came_from_other( is_assigned_to_k)) / N_k;
    falseNegGauss = sum(prob_came_from_k    (~is_assigned_to_k)) / N_k;
    
    % Compute the isolation distance and L-ratio as well
    mahalDistSq_otherSpikes = gaussMahalSq(~is_assigned_to_k, k);
    
    % Isolation distance
    mahalDistSq_sorted = sort(mahalDistSq_otherSpikes);
    
    if (N_k < length(mahalDistSq_sorted))
        isolationDist = mahalDistSq_sorted(N_k);
    else
        isolationDist = [];
    end
    
    % L-ratio
    Lratio = sum(chi2cdf(mahalDistSq_otherSpikes, nDims, 'upper')) / N_k;
    
    % Save these values
    
    metrics(k).id     = clustIDs(k) - n; % Cluster ID
    metrics(k).Lratio = Lratio;          % L-ratio
    metrics(k).IsoDis = isolationDist;   % isolation distance
    metrics(k).FP_t   = falsePosMODT;    % False positives (T-distribution)
    metrics(k).FN_t   = falseNegMODT;    % False negatives (T-distribution)
    metrics(k).FP_g   = falsePosGauss;   % False positives (Gaussian)
    metrics(k).FN_g   = falseNegGauss;   % False negatives (Gaussian)
    
end

%% Save isolation metrics

spikes   = spikesOld;
clustID1 = [metrics.id];
clustID2 = [spikes.clusters.metrics.id];
nClusts  = length(clustID1);

for iClust = nClusts:-1:1
    I = find(clustID2 == clustID1(iClust));
    if (~isempty(I))
        spikes.clusters.metrics(I).Lratio = metrics(iClust).Lratio;
        spikes.clusters.metrics(I).IsoDis = metrics(iClust).IsoDis;
        spikes.clusters.metrics(I).FP_t   = metrics(iClust).FP_t;
        spikes.clusters.metrics(I).FN_t   = metrics(iClust).FN_t;
        spikes.clusters.metrics(I).FP_g   = metrics(iClust).FP_g;
        spikes.clusters.metrics(I).FN_g   = metrics(iClust).FN_g;
    end
end

end