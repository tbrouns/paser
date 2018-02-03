function spikes = psr_sst_filter_spikes(spikes,parameters,method)

% method = 'find' : Does not immediately delete filtered spikes, but just tags them for
% possible removal later on

if (nargin < 3); method = 'find'; end

nlength = size(spikes.spiketimes,2);

if (strcmp(method,'find') && ...
        isfield(spikes,'clusters') && isfield(spikes.clusters,'rpvs'))
    spikes.clusters = rmfield(spikes.clusters,'rpvs');
end
    
if (strcmp(method,'find') || ~isfield(spikes,'delete'))
    spikes.delete = [];
end

if (strcmp(method,'delete'))
    del = false(1,nlength);
end

%% Global correlation across probes

if (parameters.filter.spikes.corr_global && isfield(spikes,'corr_global'))
    
    if (~isfield(spikes.delete,'corr_global') || length(spikes.delete.corr_global) ~= nlength)
        spikes.delete.corr_global = spikes.corr_global > parameters.filter.spikes.corr_global_thresh;
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.corr_global; end
end

%% Mean-squared error to mean waveform of cluster

if (parameters.filter.spikes.mse_cluster)
    
    if (~isfield(spikes.delete,'mse_cluster') || length(spikes.delete.mse_cluster) ~= nlength)
        spikes.mse_cluster = psr_sst_mse_cluster(spikes,parameters);
        spikes.delete.mse_cluster = spikes.mse_cluster > parameters.filter.spikes.mse_cluster_thresh;
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.mse_cluster; end
end

%% Refractory period violations

if (parameters.filter.spikes.rpvs)
    
    if (~isfield(spikes.delete,'rpvs') || length(spikes.delete.rpvs) ~= nlength)
        [spikes.clusters.rpvs,spikes.delete.rpvs] = psr_sst_filter_rpv(spikes,parameters);
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.rpvs; end
end

%% Difference in amplitude magnitude and position

if (parameters.filter.spikes.amp)
    
    if (~isfield(spikes.delete,'amp') || length(spikes.delete.amp) ~= nlength)
        spikes.delete.amp = psr_sst_filter_amp(spikes,parameters);
    end
    
    if (strcmp(method,'delete')); del = del | spikes.delete.amp; end
end

%% Delete

if (strcmp(method,'delete')); spikes = psr_sst_spike_removal(spikes,find(del),'delete'); end

end