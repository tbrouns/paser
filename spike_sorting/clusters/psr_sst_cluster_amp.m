function [chanIDs,ampAbs,ampRel,p2p] = psr_sst_cluster_amp(spikes, clustID, parameters)

threshold  = parameters.cluster.thresh * spikes.info.bgn;
signThresh = sign(mean(threshold));
spikeIDs   = ismember(spikes.assigns,clustID);
waveforms  = spikes.waveforms(spikeIDs,:,:);
waveforms  = signThresh * squeeze(mean(waveforms,1));
ampAbs     = max(waveforms); % maximum amplitude per channel
chanIDs    = find(ampAbs  > abs(threshold)); % channels that cross threshold with mean amplitude
ampRel     = ampAbs ./ abs(threshold);
p2p        = ampAbs - min(waveforms);

end