function meanSquaredErrors = psr_sst_mse_cluster(spikes,parameters)

% Convert input
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);
    
clustIDs = unique(spikes.assigns);
nClusts  = length(clustIDs);

nspikes = length(spikes.spiketimes);
meanSquaredErrors = zeros(1,nspikes,'single');

for iClust = 1:nClusts
    
    clustID = clustIDs(iClust);
        
    % Ignore smaller clusters
    spikeIDs = find(ismember(spikes.assigns,clustID));
    if (length(spikeIDs) < parameters.cluster.min_spikes); continue; end
    
    % Find average waveform
    [~,spikeIDs] = psr_sst_amp_split(spikes,clustID,parameters);
    waveMean = mean(spikes.waveforms(spikeIDs,:,:),1);

    % Calculate mean-squared error between every spike in cluster and mean waveform
    spikeIDs = find(ismember(spikes.assigns,clustID));
    waves    = spikes.waveforms(spikeIDs,:,:);
    
    % Normalize by standard deviation
    sd = mean(spikes.info.stds);
    waveMean = psr_sst_norm_waveforms(waveMean,sd);
    waves    = psr_sst_norm_waveforms(waves,   sd);
    
    % Calculate MSE
    
    MSE = bsxfun(@minus,waveMean,waves);
    MSE = mean(MSE.^2,2); % mean over all samples
    MSE =  max(MSE,[],3); % take maximum channel MSE
    
    meanSquaredErrors(spikeIDs) = single(MSE);
        
end

end