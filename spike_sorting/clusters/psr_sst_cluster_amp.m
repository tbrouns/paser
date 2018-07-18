function [chanIDs,ampAbs,ampRel,p2p,peakLocs] = psr_sst_cluster_amp(spikes, clustID, parameters)

threshold  = parameters.cluster.thresh * spikes.info.bgn;
signThresh = sign(parameters.spikes.thresh);
spikeIDs   = ismember(spikes.assigns,clustID);

waveforms = spikes.waveforms(spikeIDs,:,:);
waveforms = psr_int16_to_single(waveforms,parameters);
waveforms = signThresh * squeeze(median(waveforms,1)); % calculate median because of outliers

[ampAbs,peakLocs] = max(waveforms); % maximum amplitude over sample points per channel
chanIDs           = find(ampAbs  > abs(threshold)); % channels that cross threshold with median amplitude
ampRel            = ampAbs ./ abs(threshold);
p2p               = ampAbs - min(waveforms);

end