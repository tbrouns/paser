function spikes = psr_sst_sorting_OST(spikes,parameters)

psr_sst_ost_setpath(parameters.path.ost)

spikeWaveforms = psr_int16_to_single(spikes.waveforms(:,:),parameters);
stdEstimate = mean(spikes.info.bgn);
sortTill = 9999999;

% mainSimulatedEval;

[assigns, nrAssigned, baseSpikes, baseSpikesID] = psr_sst_ost_sortSpikesOnline(spikeWaveforms, stdEstimate, sortTill);

psr_remove_path(parameters.path.ost);
spikes.assigns = assigns;

end
