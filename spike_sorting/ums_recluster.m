function spikes = ums_recluster(spikes,clusters,metadata)

% Ignore single unit clusters

tf = strcmp({clusters.vars(:).unit},{'single'});
clusters.vars(tf) = [];

spikes.labels(tf,:) = [];

% Collect spikes from non-single units

clusterIDs = reshape([clusters.vars.id],size(clusters.vars));
spikeIDs   = ismember(spikes.assigns,clusterIDs); 

spikes.waveforms        = spikes.waveforms      (spikeIDs,:,:);
spikes.spiketimes       = spikes.spiketimes     (spikeIDs);
spikes.unwrapped_times  = spikes.unwrapped_times(spikeIDs);
spikes.trials           = spikes.trials         (spikeIDs);
spikes.amplitudes       = spikes.amplitudes     (spikeIDs);
spikes.assigns          = spikes.assigns        (spikeIDs);

spikes.info.detect.event_channel = spikes.info.detect.event_channel(spikeIDs);

spikes = ums_clustering(spikes,metadata);

spikes   = ss_aggregate(spikes);
clusters = ums_clusterfilter(spikes); %#ok

end