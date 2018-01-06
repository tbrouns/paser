function meanSquaredErrors = psr_sst_mse_cluster(spikes,parameters)

% Convert input
if (isa(spikes.waveforms,'int16'))
    spikes.waveforms = psr_single(spikes.waveforms,parameters);
end
    
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
    sd(1,1,:) = mean(spikes.info.stds);
    waveMean  = waveMean ./ repmat(sd, [size(waveMean, 1) size(waveMean, 2) 1]);
    waves     = waves    ./ repmat(sd, [size(waves,    1) size(waves,    2) 1]);
    
    % Calculate MSE
    
    MSE = bsxfun(@minus,waveMean,waves);
    MSE = mean(MSE.^2,2); % mean over all samples
    MSE = max(MSE,[],3); % take maximum channel MSE
    
    meanSquaredErrors(spikeIDs) = single(MSE);
        
end

end