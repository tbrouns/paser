function spikes = psr_sst_remove_spikes(spikes,id,method)

% Check input
nlength = length(spikes.spiketimes);

if (islogical(id))
    id = find(id); 
end

if (~isempty(id))
    id = unique(id);
    id = sort(id);
    if (id(1) < 1 || id(end) > nlength)
        disp('Error in "psr_sst_remove_spikes": spike indices are out of bounds'); 
        return;
    end
end

if (nargin < 3); method = 'delete'; end

if (strcmp(method,'delete')) % delete all input indices

    if (~isempty_field(spikes,'spikes.spiketimes'));    spikes.spiketimes   (:,id)   = []; end
    if (~isempty_field(spikes,'spikes.waveforms'));     spikes.waveforms    (id,:,:) = []; end
    if (~isempty_field(spikes,'spikes.assigns'));       spikes.assigns      (:,id)   = []; end
    if (~isempty_field(spikes,'spikes.assigns_prior')); spikes.assigns_prior(:,id)   = []; end
    if (~isempty_field(spikes,'spikes.features'));      spikes.features     (:,id)   = []; end
    if (~isempty_field(spikes,'spikes.blocks'));        spikes.blocks       (:,id)   = []; end
    
    if (~isempty_field(spikes,'spikes.delete'))
        fields = fieldnames(spikes.delete);
        for iField = 1:length(fields); spikes.delete.(fields{iField})(:,id) = []; end
    end
    
elseif (strcmp(method,'keep')) % keep only the input indices
    
    if (~isempty_field(spikes,'spikes.spiketimes'));    spikes.spiketimes    = spikes.spiketimes   (:,id);   end
    if (~isempty_field(spikes,'spikes.waveforms'));     spikes.waveforms     = spikes.waveforms    (id,:,:); end
    if (~isempty_field(spikes,'spikes.assigns'));       spikes.assigns       = spikes.assigns      (:,id);   end
    if (~isempty_field(spikes,'spikes.assigns_prior')); spikes.assigns_prior = spikes.assigns_prior(:,id);   end
    if (~isempty_field(spikes,'spikes.features'));      spikes.features      = spikes.features     (:,id);   end
    if (~isempty_field(spikes,'spikes.blocks'));        spikes.blocks        = spikes.blocks       (:,id);   end
        
    if (~isempty_field(spikes,'spikes.delete'))
        fields = fieldnames(spikes.delete);
        for iField = 1:length(fields); spikes.delete.(fields{iField}) = spikes.delete.(fields{iField})(:,id); end
    end
    
end

end