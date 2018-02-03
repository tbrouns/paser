function [cAuc,xDist,yDist] = psr_sst_cluster_stability(spikes, clustID, parameters)

% Constants

f_lambda = 5;
k_ratio  = 0.025;

% Parse inputs

spikeIDs = ismember(spikes.assigns, clustID);
spiketimes = sort(spikes.spiketimes(spikeIDs));
nspikes = length(spiketimes);
twin    = parameters.cluster.stability.win;
tlength = sum(spikes.info.dur);
frate   = nspikes / tlength;
nWin    = floor(2 * tlength / twin); % number of windows
counts  = zeros(nWin,1);
lambda  = frate * twin;

% Acquire spike counts
tstart = 1;
for i = 1:nWin
    tend = tstart + twin; % Update frequency window
    if (tend > tlength) % Check if end of signal is reached
        tstart = tlength - twin;
        tend   = tlength;
    end
    counts(i) = sum(spiketimes >= tstart & spiketimes < tend);
    tstart = tstart + round(0.5 * twin); % Update starting frequency
end

k_max  = max(counts);
if (k_max > f_lambda * lambda); k_max = f_lambda * lambda; end
k_step = round(k_ratio * k_max);
if (k_step < 1); k_step = 1; end
edges = -0.5 * k_step : k_step : k_max + 0.5 * k_step;
if (edges(end) < k_max); edges(end+1) = k_max + 0.5 * k_step; end
yDist = histcounts(counts, edges, 'Normalization', 'probability');

% Gaussian/Poisson distribution

xDist = edges(1:end-1) + 0.5 * k_step;
y = psr_sst_normpdf_stability(xDist,lambda,parameters);

% Find difference between theoretical and empirical Poisson distributions

I = yDist > y;
C = yDist;
C(I) = y(I);
cAuc = sum(C) / sum(yDist); % area under Poisson curve (relative to total area)

end