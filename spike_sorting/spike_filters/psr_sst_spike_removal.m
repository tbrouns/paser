function spikes = psr_sst_spike_removal(spikes,id,method)

% Check input
if (~isfield(spikes,'spiketimes')); disp('Error: "spikes" structure does not contain "spiketimes" field'); return; end
nlength = length(spikes.spiketimes);
if (~isempty(id))
    id = unique(id);
    id = sort(id);
    if (id(1) < 1 || id(end) > nlength); disp('Error: spike indices are out of bounds'); return; end
end

if (nargin < 3); method = 'delete'; end

if (strcmp(method,'delete')) % delete all input indices

    spikes.spiketimes(id)     = [];
    spikes.waveforms (id,:,:) = [];
            
    if (isfield(spikes,'assigns'));       spikes.assigns      (id) = []; end
    if (isfield(spikes,'assigns_prior')); spikes.assigns_prior(id) = []; end
    if (isfield(spikes,'corr_global'));   spikes.corr_global  (id) = []; end
    if (isfield(spikes,'mse_cluster'));   spikes.mse_cluster  (id) = []; end
    if (isfield(spikes,'trials'));        spikes.trials       (id) = []; end
    
    if (isfield(spikes,'delete'))
        fields = fieldnames(spikes.delete);
        for iField = 1:length(fields); spikes.delete.(fields{iField}) = []; end
    end
    
elseif (strcmp(method,'keep')) % keep only the input indices
    
    spikes.spiketimes = spikes.spiketimes(id);
    spikes.waveforms  = spikes.waveforms (id,:,:);
    
    if (isfield(spikes,'assigns'));       spikes.assigns       = spikes.assigns      (id); end
    if (isfield(spikes,'assigns_prior')); spikes.assigns_prior = spikes.assigns_prior(id); end
    if (isfield(spikes,'corr_global'));   spikes.corr_global   = spikes.corr_global  (id); end
    if (isfield(spikes,'mse_cluster'));   spikes.mse_cluster   = spikes.mse_cluster  (id); end
    if (isfield(spikes,'trials'));        spikes.trials        = spikes.trials       (id); end
        
    if (isfield(spikes,'delete'))
        fields = fieldnames(spikes.delete);
        for iField = 1:length(fields); spikes.delete.(fields{iField}) = spikes.delete.(fields{iField})(id); end
    end
    
end

end