function spikes = psr_sst_filter_chan_loc(spikes,parameters)

% Remove noise from sub-threshold channels if channel peak does not occur
% at same location as maximum spike peak

Fs = spikes.Fs;
clusterIDs = unique(spikes.assigns);
nChans = size(spikes.waveforms,3);

clustIDs_saved = [];
if (~isempty_field(spikes,'spikes.clusters.noise.id'))
    clustIDs_saved = [spikes.clusters.noise.id];
end

for iClust = fliplr(clusterIDs)
              
    % Find location of maximum waveform peak 
    chanMaxID          = psr_sst_max_amp_chan(spikes,iClust,parameters);
    [~,~,~,~,peakLocs] = psr_sst_cluster_amp (spikes,iClust,parameters);
    peakLocMax = peakLocs(chanMaxID);
    
    lowerIDs = true(1,nChans);
    lowerIDs(chanMaxID) = false;
    lowerIDs = find(lowerIDs);
    
    %% Calculate distance between peak locations 
    dLocs = zeros(nChans,1); 
    for iChan = lowerIDs
        dLocs(iChan) = abs(peakLocMax - peakLocs(iChan));
    end
    dt = 1000 * dLocs / Fs;
    I = find(dt > parameters.spikes.max_desync);
    
%     if (~isempty(I)) % TEMP
%         close all; 
%         figure; set(gcf,'position',get(0,'screensize')); 
%         
%         % Extract waveforms for this cluster
%         spikeIDs  = ismember(spikes.assigns,iClust);
%         waveforms = psr_int16_to_single(spikes.waveforms(spikeIDs,:,:),parameters);
%         waveforms = psr_sst_norm_waveforms(waveforms,spikes.info.thresh); % Normalize waveforms by spike threshold
%         waveformMedian = median(waveforms,1); % Mean across spikes
%         
%         subplot(2,1,1); hold on;
%         for iChan = 1:size(waveformMedian,3)
%             w = waveformMedian(1,:,iChan);
%             w = w(:);
%             plot(1:length(w),w);
%         end 
%         legend('Chan 1','Chan 2','Chan 3','Chan 4');
%         
%         subplot(2,1,2);
%         psr_parameters_display;
%         psr_sst_plot_waveforms(spikes, iClust, parameters);
%         
%         suplabel(num2str(dt'));
%         
%         % Save        
%         savePathFig = 'G:\Data\electrophysiology\AVP\data_figures\LOC\';
%         filename = join([spikes.session,spikes.probeID,num2str(iClust)],'-');
%         savePathFig = [savePathFig filename{1}];
%         export_fig(savePathFig);
%     end
    
    %% Save
    jClust = find(iClust == clustIDs_saved);
    if (isempty(jClust))
        if (~isempty_field(spikes,'spikes.clusters.noise')); jClust = length(spikes.clusters.noise) + 1;
        else,                                                    jClust = 1;
        end
        spikes.clusters.noise(jClust).id = iClust;
    end
    spikes.clusters.noise(jClust).loc = I;
    
end

spikes = psr_sst_white_noise(spikes,parameters);

end