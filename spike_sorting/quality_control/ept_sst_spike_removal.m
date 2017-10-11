function spikes = ept_sst_spike_removal(spikes,id,method)

% Check input
if (isempty(id)); return; end
if (~isfield(spikes,'spiketimes')); disp('Error: "spikes" structure does not contain "spiketimes" field'); return; end
nlength = length(spikes.spiketimes);
id = unique(id);
id = sort(id);
if (id(1) < 1 || id(end) > nlength); disp('Error: spike indices are out of bounds'); return; end

if (nargin < 3); method = 'delete'; end

if (strcmp(method,'delete')) % delete all input indices

    spikes.spiketimes(id)     = [];
    spikes.waveforms (id,:,:) = [];
            
    if (isfield(spikes,'assigns'));         spikes.assigns         (id) = []; end
    if (isfield(spikes,'assigns_prior'));   spikes.assigns_prior   (id) = []; end
    if (isfield(spikes,'correlations'));    spikes.correlations    (id) = []; end
    if (isfield(spikes,'removed'));         spikes.removed         (id) = []; end
    if (isfield(spikes,'trials'));          spikes.trials          (id) = []; end
    if (isfield(spikes,'unwrapped_times')); spikes.unwrapped_times (id) = []; end
    
elseif (strcmp(method,'keep')) % keep only the input indices
    
    spikes.spiketimes = spikes.spiketimes(id);
    spikes.waveforms  = spikes.waveforms (id,:,:);
    
    if (isfield(spikes,'assigns'));         spikes.assigns         = spikes.assigns         (id); end
    if (isfield(spikes,'assigns_prior'));   spikes.assigns_prior   = spikes.assigns_prior   (id); end
    if (isfield(spikes,'correlations'));    spikes.correlations    = spikes.correlations    (id); end
    if (isfield(spikes,'removed'));         spikes.removed         = spikes.removed         (id); end
    if (isfield(spikes,'trials'));          spikes.trials          = spikes.trials          (id); end
    if (isfield(spikes,'unwrapped_times')); spikes.unwrapped_times = spikes.unwrapped_times (id); end
        
elseif (strcmp(method,'array')) % don't remove, but create new array
    if (~isfield(spikes,'removed') || length(spikes.removed) ~= nlength)
        spikes.removed = zeros(size(spikes.spiketimes),'logical'); 
    end
    spikes.removed(id) = 1;
end

end