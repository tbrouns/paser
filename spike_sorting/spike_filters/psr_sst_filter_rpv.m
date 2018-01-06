function [rpvs,removed] = psr_sst_filter_rpv(spikes,parameters)

clustIDs = unique(spikes.assigns);
nClusts  = length(clustIDs);
rpvs     = zeros(nClusts,2); 
removed  = false(size(spikes.spiketimes));

for iClust = 1:nClusts
    
    clustID  = clustIDs(iClust);
    spikeIDs = find(ismember(spikes.assigns, clustID));
    nspikes  = length(spikeIDs);
    if (nspikes <= parameters.cluster.min_spikes); continue; end
    
    msErrors = spikes.mse_cluster(spikeIDs);
        
    % Find refractory period violations (RPVs)
    
    spiketimes = spikes.spiketimes(spikeIDs);
    rpvIDs     = find(diff(spiketimes) <= 0.001 * parameters.spikes.ref_period);
    num_rpvs   = length(rpvIDs);
    
    % Each RPV involves two or more spikes. We remove enough spikes to
    % resolve the RPV, where we keep the spikes that have the smallest
    % mean-squared error from cluster mean
    
    rpvIDs  = rpvIDs';
    rpvIDs  = [rpvIDs,rpvIDs+1];
    [~,del] = max(msErrors(rpvIDs),[],2);
    del = sub2ind(size(rpvIDs),(1:size(rpvIDs,1))',del);
    del = spikeIDs(rpvIDs(del));
    removed(del) = true;
    rpvs(iClust,1) = num_rpvs / nspikes;
    rpvs(iClust,2) = clustID;
       
end

end