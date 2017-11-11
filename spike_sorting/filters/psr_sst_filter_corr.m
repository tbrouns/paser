function spikes = psr_sst_filter_corr(spikes,parameters,method)

id           = spikes.correlations > parameters.spikes.artifacts_corr;
artifacts    = spikeTimes(id);

% Remove spikes based on correlation

samples = zeros(1,size(data,2));
spiketimes = spikes.spiketimes; % in seconds
spiketimes = round(spiketimes * Fs) + 1;   % in sample number
samples(spiketimes)  = 2;
samples(artifacts)   = samples(artifacts) - 1;
samples(samples < 1) = [];
id = find(samples == 2);
spikes = psr_sst_spike_removal(spikes,id,method);

end