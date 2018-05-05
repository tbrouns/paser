function spikes = psr_sst_spike_align(spikes,parameters,files)

Fs       = spikes.Fs;
clustIDs = unique(spikes.assigns);
nSpikes  = size(spikes.waveforms,1);
nPoints  = size(spikes.waveforms,2);
nChans   = size(spikes.waveforms,3);
winHalf  = round(0.5 * (nPoints - 1));

% Load raw data
dataProbe = [];
nTrials = length(files);
for iTrial = 1:nTrials
    load(files{iTrial},'ts_Spikes');
    dataProbe = [dataProbe,ts_Spikes.data];
end

waveformsNew = zeros(nSpikes,nPoints,nChans,'int16');

for iClust = clustIDs
    
    spikeIDs = ismember(spikes.assigns, iClust);
    
    spiketimes = double(spikes.spiketimes(spikeIDs));
    spiketimes = round(Fs * spiketimes)' + 1;
    
    waveforms = spikes.waveforms(spikeIDs,:,:);
    waveforms = psr_int16_to_single(waveforms,parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,spikes.info.thresh);
        
    % Find location of peak in maximum amplitude channel
    [~,~,ampRel,~] = psr_sst_cluster_amp(spikes, iClust, parameters);    
    [~, chanMaxID] = max(ampRel); 
    [~,   peakLoc] = max(mean(waveforms(:,:,chanMaxID),1),[],2);
    
    d = peakLoc - winHalf;
    winidx = d + (-winHalf:winHalf);
    
    waveformsNew(spikeIDs,:,:) = psr_sst_get_waveforms(spiketimes,dataProbe,winidx);
    
end

spikes.waveforms = waveformsNew;

end