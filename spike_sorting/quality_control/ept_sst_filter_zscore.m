function spikes = ept_sst_filter_zscore(spikes,parameters,method)


clustID = [spikes.clusters.vars.id];
nspikes = [spikes.clusters.vars.nspikes];
nclusts = length(clustID);

% Remove outliers based on z-score

for iclust = 1:nclusts
    if (nspikes(iclust) >= parameters.cluster.min_spikes)
        [x,y,z,~] = plot_distances(spikes, clustID(iclust), 1, 0);
        y         = y / sum(y);
        [~,I1]    = max(y);
        I2        = find(y(I1:end) < parameters.cluster.outlier_chi,1);
        I         = I2 + I1 - 1;
        z_thresh  = x(I);
        if (~isempty(z_thresh))
            id        = z > z_thresh;
            which     = find(spikes.assigns == clustID(iclust));
            id        = which(id);
            spikes    = ept_sst_spike_removal(spikes,id,method);
        end
    end
end

end