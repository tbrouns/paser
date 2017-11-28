function spikes = psr_sst_filter_amp(spikes,parameters,method)

if (isa(spikes.waveforms,'int16'))
    spikes.waveforms = psr_single(spikes.waveforms,parameters);
end
    
amps   = max(abs(spikes.waveforms(:,:)),[],2);
id     = find(amps > parameters.cluster.max_amplitude);
spikes = psr_sst_spike_removal(spikes,id,method);

end