function removed = psr_sst_mse_cluster(spikes,parameters)

% Convert input
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

unitIDs = unique(spikes.assigns);
nUnits  = length(unitIDs);

nspikes = length(spikes.spiketimes);
meanSquaredErrors = zeros(1,nspikes,'single');

for iUnit = 1:nUnits
    
    unitID = unitIDs(iUnit);
    
    % Ignore smaller clusters
    spikeIDs = find(ismember(spikes.assigns,unitID));
    if (length(spikeIDs) < parameters.cluster.quality.min_spikes); continue; end
    
    % Find median waveform
    [~,spikeIDs] = psr_sst_amp_split(spikes,unitID,parameters);
    waveMed = median(spikes.waveforms(spikeIDs,:,:),1);
    
    % Calculate mean-squared error between every spike in cluster and median waveform
    spikeIDs = find(ismember(spikes.assigns,unitID));
    waves    = spikes.waveforms(spikeIDs,:,:);
    
    % Normalize by the background noise
    sigma   = spikes.info.bgn;
    waveMed = psr_sst_norm_waveforms(waveMed,sigma);
    waves   = psr_sst_norm_waveforms(waves,  sigma);
        
    % Calculate MSEs
    
    meanSquaredErrors(spikeIDs) = calculateMSE(waveMed,waves);
    
end

removed = meanSquaredErrors > parameters.filter.spikes.mse.thresh;

end

function MSE = calculateMSE(wMean,w)
MSE = zeros(size(w,1),1);
if (~isempty(wMean))
    MSE = bsxfun(@minus,wMean,w);
    MSE = mean(MSE.^2,2); % mean over all samples
    MSE =  max(MSE,[],3); % take maximum channel MSE
end
MSE = single(MSE);
end