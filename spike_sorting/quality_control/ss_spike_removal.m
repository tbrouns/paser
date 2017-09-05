function spikes = ss_spike_removal(spikes,id,method)

if (nargin < 3); method = 1; end

if (method) % delete all input indices

    spikes.spiketimes     (id)     = [];
    spikes.trials         (id)     = [];
    spikes.waveforms      (id,:,:) = [];
    spikes.unwrapped_times(id)     = [];
    spikes.nspikes        (id)     = [];
    
    if (isfield(spikes,'assigns'))
        spikes.assigns(id) = [];
    end
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior(id) = [];
    end
    
else % delete everything except the input indices
    
    spikes.spiketimes      = spikes.spiketimes     (id);
    spikes.trials          = spikes.trials         (id);
    spikes.waveforms       = spikes.waveforms      (id,:,:);
    spikes.unwrapped_times = spikes.unwrapped_times(id);
    spikes.nspikes         = spikes.nspikes        (id);
    
    if (isfield(spikes,'assigns'))
        spikes.assigns = spikes.assigns(id);
    end
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior = spikes.assigns_prior(id);
    end
end

end