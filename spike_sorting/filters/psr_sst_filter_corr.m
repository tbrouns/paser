function spikes = psr_sst_filter_corr(spikes,parameters,method)

id     = find(spikes.correlations > parameters.spikes.artifacts_corr);
spikes = psr_sst_spike_removal(spikes,id,method);

end