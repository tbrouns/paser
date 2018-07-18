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
    
    % Find maximum amplitude channel
    chanMaxID = psr_sst_max_amp_chan(spikes,iClust,parameters); 
    
    % Now align each waveform to peak on the maximum amplitude channel
    waveforms = spikes.waveforms(spikeIDs,:,:);
    waveforms = psr_int16_to_single(waveforms,parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,spikes.info.thresh);
    [~,peakLocs] = max(waveforms(:,:,chanMaxID),[],2); % Peak location of each spike
    
    d = peakLocs - winHalf;
    winidx = d + (-winHalf:winHalf);
    
    spiketimes = double(spikes.spiketimes(spikeIDs));
    spiketimes = round(Fs * spiketimes)' + 1;
    waveformsNew(spikeIDs,:,:) = psr_sst_get_waveforms(spiketimes,dataProbe,winidx);
    
end

spikes.waveforms = waveformsNew;

end