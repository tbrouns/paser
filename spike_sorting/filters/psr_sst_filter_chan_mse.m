function spikes = psr_sst_filter_chan_mse(spikes,parameters)

% Remove noise from sub-threshold channels if mean-squared error of the
% normalized waveforms between the channel and spike channel is large

clusterIDs = unique(spikes.assigns);
thresh     = parameters.filter.chan.mse.thresh; % MSE threshold
nChans     = size(spikes.waveforms,3);

clustIDs_saved = [];
if (~isempty_field(spikes,'spikes.clusters.noise.id'))
    clustIDs_saved = [spikes.clusters.noise.id];
end

for iClust = fliplr(clusterIDs)
    
    % Extract cluster ID
    spikeIDs  = ismember(spikes.assigns,iClust);
    waveforms = psr_int16_to_single(spikes.waveforms(spikeIDs,:,:),parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,spikes.info.thresh); % Normalize waveforms by spike threshold
    
    % Find maximum amplitude channel
    chanMaxID = psr_sst_max_amp_chan(spikes,iClust,parameters);
    
    lowerIDs = true(1,nChans);
    lowerIDs(chanMaxID) = false;
    lowerIDs = find(lowerIDs);
    
    % Normalize waveforms by amplitude
    amplitudes = median(max(waveforms,[],2),1); % Median maximum amplitude for each channel
    waveforms = psr_sst_norm_waveforms(waveforms,amplitudes);
    waveformMedian = median(waveforms,1);
    
    %% Calculate mean-squared error (MSE) between channels
    
    MSE = zeros(nChans,1);
    for iLower = lowerIDs
        d = waveformMedian(:,:,chanMaxID) - waveformMedian(:,:,iLower);
        d = d(:);
        MSE(iLower) = mean(d.^2); % mean-squared error
    end
    I = find(MSE > thresh); % Tag channel if above threshold
    
%     if (~isempty(I)) % TEMP
%         close all; 
%         figure; set(gcf,'position',get(0,'screensize')); 
%         
%         subplot(2,1,1); hold on;
%         for iChan = 1:size(waveformMedian,3) 
%             plot(waveformMedian(1,:,iChan)); 
%         end 
%         legend('Chan 1','Chan 2','Chan 3','Chan 4');
%         title(['Chan ID max: ' num2str(chanMaxID)]);
%         
%         subplot(2,1,2);
%         psr_parameters_display;
%         psr_sst_plot_waveforms(spikes,iClust,parameters);
%         
%         suplabel(num2str(MSE'));
%         
%         % Save        
%         savePathFig = 'G:\Data\electrophysiology\AVP\data_figures\MSE\';
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
    spikes.clusters.noise(jClust).mse = I;
    
end

spikes = psr_sst_white_noise(spikes,parameters);

end