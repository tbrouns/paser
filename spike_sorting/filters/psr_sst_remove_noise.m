function spikes = psr_sst_remove_noise(spikes,parameters)

% Remove noise from sub-threshold channels if mean-squared error of the
% normalized waveforms between the channel and spike channel is large

clusterIDs = unique(spikes.assigns);
nSamples   = size(spikes.waveforms,2);
thresh     = parameters.filter.spikes.noise.thresh; % MSE threshold
precision  = 10^parameters.general.precision;
sigma      = spikes.info.bgn; 

VISUALIZE = true;
if (VISUALIZE); figure; set(gcf,'position',get(0,'screensize')); end % TEMP

for iClust = clusterIDs
    
    if (VISUALIZE)
        clf; % TEMP
        subplot(1,2,1);
        psr_sst_plot_waveforms(spikes, iClust, parameters);
    end
    
    % Extract cluster ID
    spikeIDs  = ismember(spikes.assigns,iClust);
    waveforms = psr_int16_to_single(spikes.waveforms(spikeIDs,:,:),parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,spikes.info.thresh); % Normalize waveforms by spike threshold
    
    % Find channels that exceed threshold
    amplitudes = mean(max(waveforms,[],2),1); % Mean maximum amplitude for each channel
    I = amplitudes(:)' >= 1.0;
    lowerIDs = find(~I);
    upperIDs = find( I);
    nLower = length(lowerIDs);
    nUpper = length(upperIDs);
            
    % Normalize waveforms by amplitude
    waveforms = psr_sst_norm_waveforms(waveforms,amplitudes);
    waveforms = mean(waveforms,1);
    
    % Calculate mean-squared error (MSE) between channels
    for iChan = 1:nLower
        MSE = zeros(nUpper,1);
        lowerID = lowerIDs(iChan);
        for jChan = 1:nUpper
            d = waveforms(:,:,upperIDs(jChan)) - waveforms(:,:,lowerID);
            MSE(jChan) = mean(d.^2); % mean-squared error
        end
        if (MSE > thresh) % If MSE below threshold: substitute channels with Gaussian noise
            R = normrnd(0,sigma(lowerID),sum(spikeIDs),nSamples);
            spikes.waveforms(spikeIDs,:,lowerID) = int16(precision * R);
        end
    end
    
    if (VISUALIZE)
        subplot(1,2,2); % TEMP
        psr_sst_plot_waveforms(spikes, iClust, parameters);
        export_fig(['P' num2str(parameters.probeID) '_C' num2str(iClust)]);
    end
    
end

end