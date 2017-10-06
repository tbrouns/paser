function metrics = ept_sst_sorting_MDT(spikes,parameters,weights)

% Mixture of drifting t-distribution model for sorting spikes and measuring unit isolation
% 
% This function is based on 'demo_isolation_metrics' by Kevin Shan [Copyright (c) 2016] 
% See LICENSE in MDT folder (or: https://github.com/kqshan/MoDT)
%
% Major edits by Terence Brouns (2017)

if (nargin < 3); weights = ones(size(spikes.spiketimes)); end
dims = parameters.sorting.mdt.dims;

% Remove clusters that are too small

clustIDs = unique(spikes.assigns);
nclusts  = length(clustIDs);
nspikes  = zeros(1,nclusts);
for iClust = 1:nclusts
    nspikes(iClust) = sum(spikes.assigns == clustIDs(iClust));
end

I        = nspikes <= 2 * dims;
del      = clustIDs(I);
ndel     = length(del); % clusters to remove
for iClust = 1:ndel
    id = find(spikes.assigns == del(iClust));
    spikes = ept_sst_spike_removal(spikes,id,'delete');
end
clustIDs = clustIDs(~I); % Keep remaining clusters

% Do PCA

PC = ept_pca(spikes,dims);

% Data conversion

%  New data structure with fields:
%     spk_Y           [D x N] spikes (N spikes in a D-dimensional feature space)
%     spk_t           [N x 1] spike times (ms)
%     spk_clustId     [N x 1] cluster ID this spike is assigned to
%     weighting       [N x 1] suggested weighting for training on a subset

data             = [];
data.spk_Y       = double(PC');
data.spk_t       = double(spikes.spiketimes' * 1000); % convert to ms
data.spk_clustId = double(spikes.assigns');
data.weighting   = double(weights');

% Run algorithm

% These are the parameters that we recommend in the paper
nu = parameters.sorting.mdt.nu;                               % t-distribution nu parameter (smaller = heavier tails)
q_perhour = parameters.sorting.mdt.q_perhour;                 % Drift regularization (smaller = more smoothing)
timeframe_minutes = parameters.sorting.mdt.timeframe_minutes; % Time frame duration (mostly a computational thing)

% Construct an MoDT object using these parameters
q_perframe = q_perhour * (timeframe_minutes / 60);
model = MoDT('nu',nu,'Q',q_perframe);

% Attach the data to the model
timeframe_ms = timeframe_minutes * 60e3;
model.attachData(data.spk_Y, data.spk_t,'frameDur',timeframe_ms);

% Fit the model parameters based on our spike assignments
fprintf('Fitting model based on spike assignments\n');
clustAssigned = data.spk_clustId;
model.initFromAssign(clustAssigned, 'verbose',true );

% Obtain the posterior probability that spike n came from cluster k
posterior = model.getValue('posterior');

% Let's also fit a drifting Gaussian model
fprintf('Fitting a drifting Gaussian model by setting nu=Inf\n');
gaussModel = model.copy();
gaussModel.setParams('nu',Inf);
gaussModel.initFromAssign( clustAssigned );
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
    falsePosMODT = sum(prob_came_from_other( is_assigned_to_k) ) / N_k;
    falseNegMODT = sum(prob_came_from_k    (~is_assigned_to_k) ) / N_k;
    
    % Repeat this for the Gaussian
    prob_came_from_k = gaussPosterior(:,k);
    prob_came_from_other = sum(gaussPosterior(:,otherClusts), 2);
    falsePosGauss = sum(prob_came_from_other( is_assigned_to_k)) / N_k;
    falseNegGauss = sum(prob_came_from_k    (~is_assigned_to_k)) / N_k;
    
    % Compute the isolation distance and L-ratio as well
    mahalDistSq_otherSpikes = gaussMahalSq(~is_assigned_to_k, k);
    
    % Isolation distance
    mahalDistSq_sorted = sort(mahalDistSq_otherSpikes);
    isolationDist = mahalDistSq_sorted(N_k);
    
    % L-ratio
    Lratio = sum(chi2cdf(mahalDistSq_otherSpikes, nDims, 'upper')) / N_k;
    
    % Save these values

    metrics(k).id     = clustIDs(k);   % Cluster ID
    metrics(k).Lratio = Lratio;        % L-ratio
    metrics(k).IsoDis = isolationDist; % isolation distance
    metrics(k).FP_t   = falsePosMODT;  % False positives (T-distribution)
    metrics(k).FN_t   = falseNegMODT;  % False negatives (T-distribution)
    metrics(k).FP_g   = falsePosGauss; % False positives (Gaussian)
    metrics(k).FN_g   = falseNegGauss; % False negatives (Gaussian)
    
end

end