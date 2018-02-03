function spikes = psr_sst_cluster_thresholds(spikes,parameters)

% Type definitions:
% type = 0: noise cluster
% type = 1: non-isolated single unit
% type = 2: isolated single unit

%% Hard thresholds

thresh  = parameters.cluster;
metrics = spikes.clusters.metrics;
N = size(metrics,2);

%% Check for single unit quality

I = ones(N,1);
M = {metrics.ampRel};  T = thresh.min_amp;    for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} < T; I(i) = 0; end; end % mean sub-threshold amplitude
M = {metrics.amp};     T = thresh.max_amp;    for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} > T; I(i) = 0; end; end % absolute amplitude
M = {metrics.p2p};     T = thresh.max_p2p;    for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} > T; I(i) = 0; end; end % absolute peak-to-peak
M = {metrics.sub};     T = thresh.max_sub;    for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} > T; I(i) = 0; end; end % fraction of sub-threshold spikes
M = {metrics.rpv};     T = thresh.max_rpv;    for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} > T; I(i) = 0; end; end % fraction of refractory period violations
M = {metrics.nspikes}; T = thresh.min_spikes; for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} < T; I(i) = 0; end; end % small spike number
M = {metrics.frate};   T = thresh.min_frate;  for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} < T; I(i) = 0; end; end % firing rate
M = {metrics.cAuc};    T = thresh.min_auc;    for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} < T; I(i) = 0; end; end % poisson area under curve
M = {metrics.xcLag};   T = thresh.max_xclag;  for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} > T; I(i) = 0; end; end % peak cross-correlation lag

type = I;

%% Check for isolated single unit

I = ones(N,1);
M = {metrics.Lratio};  T = thresh.max_lratio;  for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} > T; I(i) = 0; end; end % l-ratio
M = {metrics.IsoDis};  T = thresh.min_isodist; for i = 1:N; if ~isempty(M{i}) && ~isempty(T) && M{i} < T; I(i) = 0; end; end % isolation distance

type = type + type .* I;

%% Save

for i = 1:N; spikes.clusters.metrics(i).type = type(i); end

end