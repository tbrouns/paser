function spikes = ss_pca(spikes)

samples_before = round(spikes.params.Fs * spikes.params.cross_time  / 1000);
jitter_range   = samples_before - 1 + (1:round(spikes.params.max_jitter * spikes.params.Fs / 1000));

% identify which channel the event occurred on
divisor = repmat(spikes.info.detect.thresh, [size(spikes.waveforms,1) 1]);
[~, spikes.info.detect.event_channel] = max(squeeze(min(spikes.waveforms(:,jitter_range,:), [], 2))./divisor, [], 2);
spikes.info.detect.event_channel = single(spikes.info.detect.event_channel);

% save some more data that will be useful later
spikes.info.detect.align_sample = samples_before + 1;
[pca.u,pca.s,pca.v] = svd(detrend(spikes.waveforms(:,:),'constant'), 0); % SVD the data matrix
spikes.info.pca = pca;

end