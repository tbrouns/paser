function spikes = psr_sst_filter_pca(spikes,parameters,method)

clustID = [spikes.clusters.vars.id];
nClusts = length(clustID);

PC = psr_pca(spikes,parameters.cluster.pca_dims);

% Remove outliers based on PC distance

for iClust = 1:nClusts
    
    which = find(spikes.assigns == clustID(iClust));
    nspikes = length(which);
    if (nspikes > parameters.cluster.pca_dims * parameters.cluster.min_spikes)
        PC_clust = PC(which,:);
        D  = mahal(PC_clust,PC_clust);
        sd = std(D);
        id = D > sd * parameters.cluster.outlier_std;
        id = which(id);
        spikes = psr_sst_spike_removal(spikes,id,method);
    end
end

end