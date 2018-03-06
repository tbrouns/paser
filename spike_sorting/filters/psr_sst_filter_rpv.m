function removed = psr_sst_filter_rpv(spikes,parameters)

clustIDs = unique(spikes.assigns);
nClusts  = length(clustIDs);
removed  = false(size(spikes.spiketimes));

for iClust = 1:nClusts
    
    clustID  = clustIDs(iClust);
    spikeIDs = find(ismember(spikes.assigns, clustID));
    nspikes  = length(spikeIDs);
    if (nspikes <= parameters.cluster.quality.min_spikes); continue; end
        
    featureCluster = spikes.features(:,spikeIDs);
    centroid = median(featureCluster,2);
    
    % Find refractory period violations (RPVs)
    
    spiketimes = spikes.spiketimes(spikeIDs);
    rpvIDs     = find(diff(spiketimes) <= 0.001 * parameters.spikes.ref_period);
    
    % Each RPV involves two or more spikes. We remove enough spikes to
    % resolve the RPV, where we keep the spikes that have the smallest
    % mean-squared error from cluster mean
        
    nrpvs = length(rpvIDs);
    del = zeros(nrpvs,1);

    for i = 1:nrpvs
        pc = zeros(2,1);
        for j = 0:1
            pc(j+1) = sum(abs(featureCluster(:,rpvIDs(i)+j) - centroid));
        end
        [~,I] = min(pc);
        del(i) = rpvIDs(i) + I - 1;
    end

    del = unique(spikeIDs(del));
    removed(del) = true;
    
end

end