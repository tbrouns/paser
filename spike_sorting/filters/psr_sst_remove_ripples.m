function spikes = psr_sst_remove_ripples(spikes,files,parameters)

clusterIDs = unique(spikes.assigns);
nSamples   = size(spikes.waveforms,2);
nChans     = size(spikes.waveforms,3);
precision  = 10^parameters.general.precision;
sigma      = spikes.info.bgn;

minLag = spikes.Fs / parameters.filter.spikes.ripple.freq_max; 
maxLag = spikes.Fs / parameters.filter.spikes.ripple.freq_min; 
window = 2 * maxLag;
winidx = -window:window;
thresh = parameters.filter.spikes.ripple.corr_max;

VISUALIZE = false;
if (VISUALIZE) % TEMP
    close all
    psr_parameters_display;
    figure; set(gcf,'position',get(0,'screensize')); 
end 

% Load raw data
dataProbe = [];
nTrials = length(files);
for iTrial = 1:nTrials
    load(files{iTrial},'ts_Spikes');
    dataProbe = [dataProbe,ts_Spikes.data];
end

for iClust = clusterIDs
    
    if (VISUALIZE)
        clf; % TEMP
        subplot(2,2,1);
        psr_sst_plot_waveforms(spikes, iClust, parameters);
    end

    % Extract cluster ID
    spikeIDs  = ismember(spikes.assigns, iClust);
    nSpikes   = sum(spikeIDs);
    spiketimes = double(spikes.spiketimes(spikeIDs));
    spiketimes = round(spiketimes * spikes.Fs)';
    
    waveforms = psr_sst_get_waveforms(spiketimes,dataProbe,winidx);
    waveforms = mean(waveforms,1);
    waveforms = squeeze(waveforms);
    
    if (VISUALIZE)
        subplot(2,2,3);
        plot(waveforms(:));
    end
    
    % Calculate auto-correlation for individual channels
    % Large peaks at non-zero lag indicative of noise ripples
    autoCorr = [];
    for iChan = nChans:-1:1
        autoCorr(:,iChan) = xcorr(waveforms(:,iChan),maxLag,'coeff');
    end
    
    % Ignore points around zero lag
    d = maxLag - minLag;
    autoCorr = autoCorr([1:d,end-d+1:end],:);
    
    % Threshold
    autoCorr = max(autoCorr);
    I = find(autoCorr > thresh);
    
    if (any(I)) % Replace with white noise
        for iChan = I
            R = normrnd(0,sigma(iChan),sum(spikeIDs),nSamples);
            spikes.waveforms(spikeIDs,:,iChan) = int16(precision * R);
        end
    end
    
    if (VISUALIZE)
        subplot(2,2,2); % TEMP
        psr_sst_plot_waveforms(spikes, iClust, parameters);
        export_fig(['P' num2str(parameters.probeID) '_C' num2str(iClust)]);
    end
    
end

end