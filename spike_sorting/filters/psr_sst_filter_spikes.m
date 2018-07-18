function spikes = psr_sst_filter_spikes(spikes,parameters,method)

% method = 'find' : Does not immediately delete filtered spikes, but just tags them for
% possible removal later on

if (nargin < 3); method = 'find'; end

nlength = size(spikes.spiketimes,2);

if (strcmp(method,'find') && ~isempty_field(spikes,'spikes.clusters.rpvs'))
    spikes.clusters = rmfield(spikes.clusters,'rpvs');
end

if (~psr_isfield(spikes,'spikes.delete')); spikes.delete = []; end

switch method
    
    case 'find'
        spikes.delete = psr_remove_field(spikes.delete,'art');
        spikes.delete = psr_remove_field(spikes.delete,'mse_cluster');
        spikes.delete = psr_remove_field(spikes.delete,'rpvs');
        spikes.delete = psr_remove_field(spikes.delete,'amp');
        
    case 'delete'
        del = false(1,nlength);
end

%% Spikes on raw signal artifacts

if (parameters.filter.spikes.art.run)
   
    if (~isfield(spikes.delete,'art') || length(spikes.delete.art) ~= nlength)
        spikes.delete.art = psr_sst_filter_art(spikes);
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.art; end
end

%% Mean-squared error to mean waveform of cluster

if (parameters.filter.spikes.mse.run)
    
    if (~isfield(spikes.delete,'mse_cluster') || length(spikes.delete.mse_cluster) ~= nlength)
        spikes.delete.mse_cluster = psr_sst_mse_cluster(spikes,parameters);
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.mse_cluster; end
end

%% Refractory period violations

if (parameters.filter.spikes.rpv.run)
    
    if (~isfield(spikes.delete,'rpvs') || length(spikes.delete.rpvs) ~= nlength)
        spikes.delete.rpvs = psr_sst_filter_rpv(spikes,parameters);
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.rpvs; end
end

%% Difference in amplitude position

if (parameters.filter.spikes.amp.run)
    
    if (~isfield(spikes.delete,'amp') || length(spikes.delete.amp) ~= nlength)
        spikes.delete.amp = psr_sst_filter_amp(spikes,parameters);
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.amp; end
end

%% Delete

if (strcmp(method,'delete')); spikes = psr_sst_remove_spikes(spikes,find(del),'delete'); end

end