function spikes = psr_sst_filter_spikes(spikes,parameters,method)

% Does not immediately delete filtered spikes, but just tags them for
% possible removal later on

if (nargin < 3); method = 'find'; end

if (strcmp(method,'find')) % find which spikes to delete
    
    spikes.delete = []; % reset
    
    % Global correlation across probes
    if (parameters.filter.spikes.corr_global)
        spikes.delete.corr_global = spikes.corr_global > parameters.filter.spikes.corr_thresh_global;
    end

    % Mean-squared error to mean waveform of cluster
    if (parameters.filter.spikes.mse_cluster)
        spikes.mse_cluster = psr_sst_mse_cluster(spikes,parameters);
        spikes.delete.mse_cluster = spikes.mse_cluster > parameters.filter.spikes.mse_thresh;
    end

    % Refractory period violations
    if (parameters.filter.spikes.rpvs) 
        [spikes.clusters.rpvs,spikes.delete.rpvs] = psr_sst_filter_rpv(spikes,parameters);
    end

    % Difference in amplitude magnitude and position
    if (parameters.filter.spikes.amp)
        spikes.delete.amp = psr_sst_filter_amp(spikes,parameters); 
    end

elseif (strcmp(method,'delete')) % delete the spikes
    if (isfield(spikes,'delete'))
        
        del = false(size(spikes.spiketimes));
        
        if (parameters.filter.spikes.corr_global && isfield(spikes.delete,'corr_global'))
            del = del | spikes.delete.corr_global;
        end
        
        if (parameters.filter.spikes.mse_cluster && isfield(spikes.delete,'mse_cluster'))
            del = del | spikes.delete.mse_cluster;
        end
        
        if (parameters.filter.spikes.rpvs && isfield(spikes.delete,'rpvs'))
            del = del | spikes.delete.rpvs;
        end
        
        if (parameters.filter.spikes.amp && isfield(spikes.delete,'amp'))
            del = del | spikes.delete.corr_global;
        end
        
        spikes = psr_sst_spike_removal(spikes,find(del),'delete');
    end

end