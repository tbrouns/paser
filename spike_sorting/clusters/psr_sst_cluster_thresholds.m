function spikes = psr_sst_cluster_thresholds(spikes,parameters)

% Type definitions:
% type = 0: noise cluster
% type = 1: multi-unit or semi-single unit
% type = 2: non-isolated single unit
% type = 3: isolated single unit

if (isempty_field(spikes,'spikes.clusters.metrics')); return; end

nClusts = size(spikes.clusters.metrics,2);

%% Initial basic amplitude check

I = ones(nClusts,1);
if psr_isfield(spikes,'spikes.clusters.metrics.ampRel'); I = checkThreshold(I, {spikes.clusters.metrics.ampRel}, parameters.cluster.quality.min_amp, true);  end % mean sub-threshold amplitude
if psr_isfield(spikes,'spikes.clusters.metrics.amp');    I = checkThreshold(I, {spikes.clusters.metrics.amp},    parameters.cluster.quality.max_amp, false); end % absolute amplitude
if psr_isfield(spikes,'spikes.clusters.metrics.p2p');    I = checkThreshold(I, {spikes.clusters.metrics.p2p},    parameters.cluster.quality.max_p2p, false); end % absolute peak-to-peak
quality = I;

%% Check for single unit quality

I = ones(nClusts,1);
if psr_isfield(spikes,'spikes.clusters.metrics.nspikes'); I = checkThreshold(I, {spikes.clusters.metrics.nspikes}, parameters.cluster.quality.min_spikes, true);  end % small spike number
if psr_isfield(spikes,'spikes.clusters.metrics.sub');     I = checkThreshold(I, {spikes.clusters.metrics.sub},     parameters.cluster.quality.max_sub,    false); end % fraction of sub-threshold spikes
if psr_isfield(spikes,'spikes.clusters.metrics.rpv');     I = checkThreshold(I, {spikes.clusters.metrics.rpv},     parameters.cluster.quality.max_rpv,    false); end % fraction of refractory period violations
if psr_isfield(spikes,'spikes.clusters.metrics.mse');     I = checkThreshold(I, {spikes.clusters.metrics.mse},     parameters.cluster.quality.max_mse,    false); end % stability measure
if psr_isfield(spikes,'spikes.clusters.metrics.xcLag');   I = checkThreshold(I, {spikes.clusters.metrics.xcLag},   parameters.cluster.quality.max_xclag,  false); end % peak cross-correlation lag
quality = quality + (quality > 0) .* I;

%% Check for isolated single unit

I = ones(nClusts,1);
if psr_isfield(spikes,'spikes.clusters.metrics.F1_t'); I = checkThreshold(I, {spikes.clusters.metrics.F1_t}, parameters.cluster.quality.min_f1, true); end % false positive rate
quality = quality + (quality > 0) .* I;

%% Save

for i = 1:nClusts; spikes.clusters.metrics(i).quality = quality(i); end

end

function I = checkThreshold(I,metric,threshold,tf)

nClusts = length(I);
for i = 1:nClusts
    if ~isempty(metric{i})
        if ~isempty(threshold)
            if (tf); if metric{i} < threshold; I(i) = 0; end
            else,    if metric{i} > threshold; I(i) = 0; end
            end
        end
    else
        I(i) = 0;
    end
end
end
