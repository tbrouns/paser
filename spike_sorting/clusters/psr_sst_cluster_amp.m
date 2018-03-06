function [chanIDs,ampAbs,ampRel,p2p] = psr_sst_cluster_amp(spikes, clustID, parameters)

threshold  = parameters.cluster.thresh * spikes.info.bgn;
signThresh = sign(mean(threshold));
spikeIDs   = ismember(spikes.assigns,clustID);
waveforms  = spikes.waveforms(spikeIDs,:,:);
waveforms  = signThresh * squeeze(mean(waveforms,1));
amplitudes = max(waveforms); % maximum amplitude per channel
chanIDs    = find(amplitudes  > abs(threshold)); % channels that cross threshold with mean amplitude
ampRel     =  max(amplitudes ./ abs(threshold));
ampAbs     =  max(amplitudes); % maximum mean channel amplitude
p2p        = ampAbs - min(waveforms(:));

end