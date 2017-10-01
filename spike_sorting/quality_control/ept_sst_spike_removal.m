function spikes = ept_sst_spike_removal(spikes,id,method)

if (nargin < 3); method = 'delete'; end

if (strcmp(method,'delete')) % delete all input indices

    spikes.spiketimes     (id)     = [];
    spikes.waveforms      (id,:,:) = [];
            
    if (isfield(spikes,'assigns'))
        spikes.assigns(id) = [];
    end
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior(id) = [];
    end
    
elseif (strcmp(method,'keep')) % keep only the input indices
    
    spikes.spiketimes      = spikes.spiketimes     (id);
    spikes.waveforms       = spikes.waveforms      (id,:,:);
    
    if (isfield(spikes,'assigns'))
        spikes.assigns = spikes.assigns(id);
    end
    
    if (isfield(spikes,'assigns_prior'))
        spikes.assigns_prior = spikes.assigns_prior(id);
    end
    
elseif (strcmp(method,'array')) % don't remove, but create new array
    if (isfield(spikes,'removed')); spikes.removed = [spikes.removed,single(id)]; 
    else                            spikes.removed = single(id);
    end
end

end