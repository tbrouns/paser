function latencies = psr_ft_latency(spikesFT)

% PSR_FT_LATENCY - Calculates latency of first spike after stimulus
%
% Syntax:  latencies = psr_ft_latency(spikesFT)
%
% Inputs:
%    spikesFT - A FieldTrip spike structure
%
% Outputs:
%    latencies - Cell array with one cell for each unit, containing a
%                vector with the latencies [sec] for all trials for that
%                unit

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

nUnits = length(spikesFT.time);
latencies = cell(nUnits,1);

for iUnit = 1:nUnits
    trials = spikesFT.trial{iUnit};
    times  = spikesFT.time {iUnit};
    nTrials = size(spikesFT.trialtime,1);
    latencies{iUnit} = NaN(nTrials,1);
    for iTrial = 1:nTrials
        timesTrial = times(trials == iTrial);
        I = find(timesTrial >= 0,1,'first');
        if ~isempty(I); latencies{iUnit}(iTrial) = timesTrial(I); end
    end
end

end