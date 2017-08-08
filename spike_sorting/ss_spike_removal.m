function spikes = ss_spike_removal(spikes,id)

    spikes.assigns        (id)     = [];
    spikes.spiketimes     (id)     = [];
    spikes.trials         (id)     = [];
    spikes.waveforms      (id,:,:) = [];
    spikes.unwrapped_times(id)     = [];
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior  (id) = [];
    end
    
end