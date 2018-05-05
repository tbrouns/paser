function latencies = psr_ft_latency(spikesFT)

nClusts = length(spikesFT.time);
latencies = cell(nClusts,1);

for iClust = 1:nClusts
    trials = spikesFT.trial{iClust};
    times  = spikesFT.time {iClust};
    nTrials = size(spikesFT.trialtime,1);
    latencies{iClust} = NaN(nTrials,1);
    for iTrial = 1:nTrials
        timesTrial = times(trials == iTrial);
        I = find(timesTrial >= 0,1,'first');
        if ~isempty(I); latencies{iClust}(iTrial) = timesTrial(I); end
    end
end

end