function spikes = psr_sst_sorting_OST(spikes,parameters)

addpath(parameters.path.ost);
setpath;

spikeWaveforms = psr_int16_to_single(spikes.waveforms(:,:),parameters);
stdEstimate = mean(spikes.info.bgn);
sortTill = 9999999;

% mainSimulatedEval;

[assigns, nrAssigned, baseSpikes, baseSpikesID] = sortSpikesOnline(spikeWaveforms, stdEstimate, sortTill);

psr_remove_path(parameters.path.ost);
spikes.assigns = assigns;

end
