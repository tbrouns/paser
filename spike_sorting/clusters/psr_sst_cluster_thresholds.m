function spikes = psr_sst_cluster_thresholds(spikes,parameters)

% Type definitions:
% type = 0: noise cluster
% type = 1: multi-unit or semi-single unit
% type = 2: non-isolated single unit
% type = 3: isolated single unit

%% Hard thresholds

nClusts = size(spikes.clusters.metrics,2);

%% Initial basic amplitude check

I = ones(nClusts,1);
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.ampRel'); I = checkThreshold({spikes.clusters.metrics.ampRel}, parameters.cluster.quality.min_amp, nClusts, I,  true); end % mean sub-threshold amplitude
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.amp');    I = checkThreshold({spikes.clusters.metrics.amp},    parameters.cluster.quality.max_amp, nClusts, I, false); end % absolute amplitude
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.p2p');    I = checkThreshold({spikes.clusters.metrics.p2p},    parameters.cluster.quality.max_p2p, nClusts, I, false); end % absolute peak-to-peak
quality = I;

%% Check for single unit quality

I = ones(nClusts,1);
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.nspikes'); I = checkThreshold({spikes.clusters.metrics.nspikes}, parameters.cluster.quality.min_spikes, nClusts, I,  true); end % small spike number
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.sub');     I = checkThreshold({spikes.clusters.metrics.sub},     parameters.cluster.quality.max_sub,    nClusts, I, false); end % fraction of sub-threshold spikes
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.rpv');     I = checkThreshold({spikes.clusters.metrics.rpv},     parameters.cluster.quality.max_rpv,    nClusts, I, false); end % fraction of refractory period violations
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.cAuc');    I = checkThreshold({spikes.clusters.metrics.cAuc},    parameters.cluster.quality.min_auc,    nClusts, I,  true); end % stability area under curve
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.xcLag');   I = checkThreshold({spikes.clusters.metrics.xcLag},   parameters.cluster.quality.max_xclag,  nClusts, I, false); end % peak cross-correlation lag
quality = quality + (quality > 0) .* I;

%% Check for isolated single unit

I = ones(nClusts,1);
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.FP_t'); I = checkThreshold({spikes.clusters.metrics.FP_t}, parameters.cluster.quality.max_fp, nClusts, I, false); end % false positive rate
if ~psr_isempty_field(spikes,'spikes.clusters.metrics.FN_t'); I = checkThreshold({spikes.clusters.metrics.FN_t}, parameters.cluster.quality.max_fn, nClusts, I, false); end % false negative rate
quality = quality + (quality > 0) .* I;

%% Save

for i = 1:nClusts; spikes.clusters.metrics(i).quality = quality(i); end

end

function I = checkThreshold(metric,threshold,nclusts,I,tf)

for i = 1:nclusts
    if ~isempty(metric{i}) && ~isempty(threshold)
        if (tf); if metric{i} < threshold; I(i) = 0; end 
        else,    if metric{i} > threshold; I(i) = 0; end 
        end
    end
end

end
