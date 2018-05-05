function ids = psr_get_spike_ids(spikePoints,dataPoints)

% spikePoints: location of spikes, given in sample number
% dataPoints: data points for which we want to look for spikes

[spikePoints,Ia,Ic] = unique(spikePoints); % Deal with spikes at same location
N = max(spikePoints);
dataPoints = dataPoints(dataPoints <= N);
X = zeros(N,1); % Binary vector of spike points
X(spikePoints) = 1;

dataPoints    = dataPoints(X(dataPoints) == 1);
X(dataPoints) = X(dataPoints) + 1;
dataPoints    = X(X > 0) == 2;
dataPoints    = Ia(dataPoints);

ids = false(size(Ic))';
ids(dataPoints) = true;
ids = ids(Ic);

end