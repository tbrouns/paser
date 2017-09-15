function spikes = ss_spikefilter_amp(spikes)

% Remove outliers based on absolute amplitude

nspikes   = size(spikes.waveforms,1);
nsamples  = size(spikes.waveforms,2);
nchans    = size(spikes.waveforms,3);
threshold = spikes.params.outlier.abs;

waves  = reshape(spikes.waveforms,nspikes,nsamples*nchans);
id     = max(sign(threshold) * waves,[],2) >= sign(threshold) * threshold;
spikes = ss_spike_removal(spikes,~id);

end