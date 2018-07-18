function spikes = psr_sst_cluster_remove(spikes)

%% Remove spikes that are part of noise clusters

if ~isempty_field(spikes,'spikes.clusters.metrics.quality')
    
    clusterIDs     = {spikes.clusters.metrics.id};
    clusterQuality = {spikes.clusters.metrics.quality};
    nClusts = length(clusterIDs);
    removed = false(size(spikes.spiketimes));
    for iClust = 1:nClusts
        if (clusterQuality{iClust} == 0)
            spikeIDs = ismember(spikes.assigns,clusterIDs{iClust});
            removed(spikeIDs) = true;
        end
    end
    
    spikes = psr_sst_remove_spikes(spikes,find(removed),'delete');
    
end

end