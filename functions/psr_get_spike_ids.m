function ids = psr_get_spike_ids(spikePoints,dataPoints)

% PSR_GET_SPIKE_IDS - Get spike IDs from given data points
%
% Syntax:  ids = psr_get_spike_ids(spikePoints,dataPoints)
%
% Inputs:
%    spikePoints - Location of spikes, given in sample number of the signal
%    dataPoints  - Data sample points for which we want to look for spikes
%
% Outputs:
%    ids - Spike IDs referencing spikes in the "spikes" structure

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

ids = [];
if (~isempty(spikePoints))
    [spikePoints,~,Ic] = unique(spikePoints); % Deal with spikes at same location
    N = max(spikePoints);
    dataPoints = dataPoints(dataPoints <= N);
    X = zeros(N,1); % Binary vector of spike points
    X(spikePoints) = 1;

    dataPoints    = dataPoints(X(dataPoints) == 1); % Data points that have a spike
    X(dataPoints) = X(dataPoints) + 1;
    dataPoints    = X(X > 0) == 2; 

    ids = false(size(Ic))';
    ids(dataPoints) = true;
    ids = ids(Ic);
end

end