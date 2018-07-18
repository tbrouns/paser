function spikes = psr_sst_filter_chan_rip(spikes,parameters,files)

Fs       = spikes.Fs;
clustIDs = unique(spikes.assigns);
nChans   = size(spikes.waveforms,3);

% The max lag shouldn't be larger than the refractory period because then
% we would erroneously remove physiological spikes that are close together. We only
% want to remove noise oscillations.

maxLag = spikes.Fs / parameters.filter.chan.rip.fmin; 
window = 2 * maxLag;
winidx = -window:window;
thresh = parameters.filter.chan.rip.thresh;

clustIDs_saved = [];
if (~isempty_field(spikes,'spikes.clusters.noise.id'))
    clustIDs_saved = [spikes.clusters.noise.id];
end

% Load raw data
nBlocks = length(files);
dataProbe = cell(1,nBlocks);
for iBlock = 1:nBlocks
    load(files{iBlock},'ts_Spikes');
    dataProbe{iBlock} = ts_Spikes.data;
end
dataProbe = cell2mat(dataProbe);
    
for iClust = clustIDs
        
    % Extract cluster ID
    spikeIDs = ismember(spikes.assigns, iClust);
    
    spiketimes = double(spikes.spiketimes(spikeIDs));
    spiketimes = round(Fs * spiketimes)' + 1;
    
    waveforms = psr_sst_get_waveforms(spiketimes,dataProbe,winidx);
    waveformMedian = median(waveforms,1);
    waveformMedian = squeeze(waveformMedian);
        
    % Calculate auto-correlation for individual channels
    % Large peaks at non-zero lag indicative of noise ripples

    for iChan = nChans:-1:1
        w = waveformMedian(:,iChan);
        Y = fft(w);
        L = length(w);
        P2 = abs(Y/L);
        P1 = P2(1:round(L/2+1));
        P1(2:end-1) = 2*P1(2:end-1);
        Y_fft(:,iChan) = P1;
    end
    
    % Threshold
    maxPower = max(Y_fft ./ sum(Y_fft));
    I = find(maxPower > thresh);
     
    %% Save
    jClust = find(iClust == clustIDs_saved);
    if (isempty(jClust))
        if (~isempty_field(spikes,'spikes.clusters.noise')); jClust = length(spikes.clusters.noise) + 1;
        else,                                                    jClust = 1;
        end
        spikes.clusters.noise(jClust).id = iClust;
    end
    spikes.clusters.noise(jClust).rip = I;
    
%     %% VISUALIZE
%     
%     close all; 
%     figure; set(gcf,'position',get(0,'screensize'));
%     
%     f = Fs*(0:size(Y_fft,1)-1);
%     for iChan = 1:nChans 
%         axs(iChan) = subplot(4,2,iChan);
%         plot(f,Y_fft(:,iChan) / sum(Y_fft(:,iChan)));
%         xlabel('Frequency (Hz)')
%         title(num2str(maxPower(iChan)));
%     end 
%     linkaxes(axs);
%     
%     subplot(4,2,[5 6]); hold on;
%     for iChan = 1:size(waveformMedian,2) 
%         plot(waveformMedian(:,iChan)); 
%     end 
%     legend('Chan 1','Chan 2','Chan 3','Chan 4');
%     
%     subplot(4,2,[7 8]);
%     psr_parameters_display;
%     psr_sst_plot_waveforms(spikes,iClust,parameters);
% 
%     % Save        
%     savePathFig = 'G:\Data\electrophysiology\AVP\data_figures\RIP\';
%     filename = join([spikes.session,spikes.probeID,num2str(iClust)],'-');
%     savePathFig = [savePathFig filename{1}];
%     export_fig(savePathFig);
    
end

spikes = psr_sst_white_noise(spikes,parameters);

end