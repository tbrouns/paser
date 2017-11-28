function spikes = psr_sst_noise_removal(spikes,parameters)

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);
nSamples   = size(spikes.waveforms,2);
nChans     = size(spikes.waveforms,3);
thresh     = parameters.spikes.corr_thresh;
precision  = 10^parameters.general.precision;

for iClust = 1:nClust
    
    % Extract cluster ID
    
    clusterID = clusterIDs(iClust);
    which     = spikes.assigns == clusterID;
    waves     = psr_single(spikes.waveforms(which,:,:),parameters);
    wMax      = mean(max(waves,[],2),1);
    [~,iMax]  = max(wMax,[],3);
    
    wNorm     = max(waves,[],2);
    wavesNorm = bsxfun(@rdivide,waves,wNorm); % Normalize
    
    channels  = zeros(nChans,1);
    
    for iChan = 1:nChans
        
        corrPair = zeros(nSamples,1);
        
        for iSample = 1:nSamples
            C = corr([wavesNorm(:,iSample,iMax),wavesNorm(:,iSample,iChan)]);
            corrPair(iSample) = C(1,2);
        end
        
        C = mean(corrPair);
        channels(iChan) = C > thresh;
        
    end
      
    noiseWaves = waves(:,:,~channels); % Channels that do not correlate with max amplitude channels
    noiseWaves = mean(noiseWaves,3);
        
    % Subtract noise from all waveforms
    
    waves = bsxfun(@minus,waves,noiseWaves);  
    spikes.waveforms(which,:,:) = int16(precision * waves);
   
end

end