function spikes = ss_spikefilter_zscore(spikes,clusters)

clustID = cell2mat({clusters.vars.id});
flags   = cell2mat({clusters.vars.flag});
clustID = clustID(flags);

nclusts = length(clustID);

% Remove outliers based on z-score

for iclust = 1:nclusts
    [x,y,z,~] = plot_distances(spikes, clustID(iclust), 1, 0);
    y         = y / sum(y);
    [~,I1]    = max(y);
    I2        = find(y(I1:end) < spikes.params.outlier.chi,1);
    I         = I2 + I1 - 1;
    z_thresh  = x(I);
    if (~isempty(z_thresh))
        id        = z > z_thresh;
        which     = find(spikes.assigns == clustID(iclust));
        id        = which(id);
        spikes    = ss_spike_removal(spikes,id);
    end
end

end