function spikes = ss_spike_removal(spikes,id,method)

if (nargin < 3); method = 1; end

if (method == 1) % delete all input indices

    spikes.spiketimes     (id)     = [];
    spikes.trials         (id)     = [];
    spikes.waveforms      (id,:,:) = [];
    spikes.unwrapped_times(id)     = [];
    spikes.nlength        (id)     = [];
    
    if (isfield(spikes,'assigns'))
        spikes.assigns(id) = [];
    end
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior(id) = [];
    end
    
elseif (method == 0) % delete everything except the input indices
    
    spikes.spiketimes      = spikes.spiketimes     (id);
    spikes.trials          = spikes.trials         (id);
    spikes.waveforms       = spikes.waveforms      (id,:,:);
    spikes.unwrapped_times = spikes.unwrapped_times(id);
    spikes.nlength         = spikes.nlength        (id);
    
    if (isfield(spikes,'assigns'))
        spikes.assigns = spikes.assigns(id);
    end
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior = spikes.assigns_prior(id);
    end
else % don't remove, but create new array
    if (isfield(spikes,'removed')); spikes.removed = [spikes.removed;single(id)]; 
    else                            spikes.removed = single(id);
    end
end

end