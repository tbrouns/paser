function spikes = psr_sst_cluster_quality(spikes,parameters)

    parameters.filter.spikes.rpv.run = false; % Don't do this yet
    
    spikes = psr_sst_filter_spikes     (spikes,parameters,'find');   % Find noise spikes ...
    spikes = psr_sst_filter_spikes     (spikes,parameters,'delete'); % ... and then delete them
    spikes = psr_sst_cluster_features  (spikes,parameters);          % Calculate cluster metrics
    spikes = psr_sst_cluster_thresholds(spikes,parameters);          % Classify cluster quality
    
end