function amplitudes = psr_sst_amp(spikes,clustID,parameters)

waveforms  = spikes.waveforms(spikes.assigns == clustID,:,:);
waveforms  = psr_int16_to_single(waveforms,parameters);
waveforms  = psr_sst_norm_waveforms(waveforms,mean(spikes.info.thresh));
amplitudes = max(waveforms(:,:),[],2);

end