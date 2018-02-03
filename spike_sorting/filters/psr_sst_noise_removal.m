function spikes = psr_sst_noise_removal(spikes,parameters)

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);
nSamples   = size(spikes.waveforms,2);
nChans     = size(spikes.waveforms,3);
thresh     = parameters.filter.spikes.mse_chan_thresh;
precision  = 10^parameters.general.precision;
sigma      = mean(spikes.info.stds); % Average over all trials

for iClust = 1:nClust
    
    % Extract cluster ID
    clustID   = clusterIDs(iClust);
    spikeIDs  = spikes.assigns == clustID;
    waveforms = psr_int16_to_single(spikes.waveforms(spikeIDs,:,:),parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,mean(spikes.info.thresh)); % Normalize waveforms by threshold
    
    % Find maximum amplitude channel
    amplitudes = mean(max(waveforms,[],2),1); % Mean maximum amplitude for each channel
    [~,iChanMax] = max(amplitudes,[],3); % Maximum amplitude over all channels
    
    % Normalize waveforms by amplitude
    waveforms = psr_sst_norm_waveforms(waveforms,amplitudes);
        
    waveforms   = mean(waveforms,1);
    waveformMax = waveforms(:,:,iChanMax);
    
    % Calculate mean-squared error (MSE) between max amplitude channel and other channels
    for iChan = 1:nChans

        d = waveformMax - waveforms(:,:,iChan); % error
        MSE = mean(d.^2); % mean-squared error
                
        % If MSE below threshold: substitute channels with Gaussian noise
        if (MSE > thresh)
            R = normrnd(0,sigma(iChan),sum(spikeIDs),nSamples);
            spikes.waveforms(spikeIDs,:,iChan) = int16(precision * R);
        end
        
    end
    
end

end