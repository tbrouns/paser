function amplitudes = psr_sst_amp(spikes,clustID,parameters)

th(1,1,:)  = mean(spikes.info.thresh);
waves      = spikes.waveforms(spikes.assigns == clustID,:,:); clear spikes;
if (isa(waves,'int16')); waves = psr_single(waves,parameters); end
waves      = waves ./ repmat(th, [size(waves,1) size(waves,2) 1]);
amplitudes = max(waves(:,:),[],2);

end