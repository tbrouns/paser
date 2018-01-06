function [SpikeTimes] = psr_get_spike_times(spikes)

clusterIDs = [spikes.clusters.metrics.id];
numclusts  = length(clusterIDs);
SpikeTimes =  cell(1,numclusts);

for iClus = 1:numclusts
    id = (spikes.assigns == clusterIDs(iClus));
    SpikeTimes{iClus} = spikes.spiketimes(id) - spikes.info.trialonset;
end

end