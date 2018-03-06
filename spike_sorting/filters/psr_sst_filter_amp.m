function removed = psr_sst_filter_amp(spikes,parameters)

% Maximum amplitude should be found in same channels as in mean waveform

spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

Fs               = spikes.Fs;
window           = round(0.5 * Fs * (parameters.spikes.max_desync / 1000));
waveforms        = spikes.waveforms;
waveforms        = psr_sst_norm_waveforms(waveforms,spikes.info.thresh);
[waveforms,Imax] = max(waveforms,[],2);
[~,maxChanIds]   = max(waveforms,[],3);

Imax = squeeze(Imax);
waveforms = squeeze(waveforms);
% I    = sub2ind(size(Imax),1:size(Imax,1),chanIDs');
% Imax = Imax(I);

clustIDs = unique(spikes.assigns);
nClusts  = length(clustIDs);
removed  = false(size(spikes.spiketimes));

for iClust = 1:nClusts
        
    % Check if max amplitude occurs in channel that exceed threshold on
    % average and if amplitude occurs at same position
    
    clustID = clustIDs(iClust);
    spikeIDs_1 = find(ismember(spikes.assigns,clustID));
    if (length(spikeIDs_1) < parameters.cluster.quality.min_spikes); continue; end
    [~,spikeIDs_2] = psr_sst_amp_split(spikes,clustID,parameters);
    ImaxCluster = Imax(spikeIDs_1,:); % Sample number of amplitude
    maxChanIdClust = maxChanIds(spikeIDs_1);
    
    chanAmps = mean(waveforms(spikeIDs_2,:));
    maxChanAmp = max(chanAmps);
    
    delChan = false(size(spikeIDs_1,2),1);
    
    if (maxChanAmp >= 1) % Should exceed threshold
        
        pksChanIdClust = find(chanAmps > maxChanAmp - parameters.filter.spikes.amp.offset);
        nChan = length(pksChanIdClust);
        
        if (nChan > 0)

            N1 = zeros(size(maxChanIdClust));
            N2 = zeros(size(maxChanIdClust));

            for iChan = 1:nChan

                % Does spike maximum occur on correct channel? 
                pksChanId = pksChanIdClust(iChan);
                I = (maxChanIdClust == pksChanId);
                N1 = N1 + I; % Correct channel if greater than zero

                % Does spike maximum occur at correct sample number?
                ImaxMode = mode(Imax(spikeIDs_2,pksChanId));
                ImaxClusChan = ImaxCluster(:,pksChanId);
                I = ImaxClusChan > ImaxMode + window | ImaxClusChan < ImaxMode - window;
                N2 = N2 + I; % Wrong position if greater than zero

            end

            delChan = (N1 == 0 | N2 ~= 0);
            
        end
    end
    
    removed(spikeIDs_1) = delChan;
    
end

end