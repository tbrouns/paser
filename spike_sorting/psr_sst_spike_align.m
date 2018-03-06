function spikes = psr_sst_spike_align(spikes,parameters)

if (parameters.filter.spikes.ripple.process)
    parameters = psr_load_parameters(parameters); % Re-load parameters to get original window size
end

precision = 10^parameters.general.precision;
win       = round(spikes.Fs * parameters.spikes.window_size / 1000);
winHalf   = round(0.5 * (win - 1));
clustIDs  = unique(spikes.assigns);
nSpikes   = size(spikes.waveforms,1);
nPoints   = size(spikes.waveforms,2);
nChans    = size(spikes.waveforms,3);

waveformsNew = NaN(nSpikes,win,nChans);

for iClust = clustIDs
    
    spikeIDs  = ismember(spikes.assigns,iClust);
    waveforms = spikes.waveforms(spikeIDs,:,:);
    waveforms = psr_int16_to_single(waveforms,parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,spikes.info.thresh);
        
    [~,locs] = max(mean(waveforms,1),[],2);
%     chanIDs = find(pks > abs(1 / parameters.spikes.thresh))'; % threshold of 1 s.d
    
    for iChan = 1:nChans
        
        loc = locs(iChan); % peak location

        d  = 0;
        i1 = loc - winHalf;
        i2 = loc + winHalf;

        if (i1 < 1);       d = i1 - 1;       i1 = 1;       end
        if (i2 > nPoints); d = i2 - nPoints; i2 = nPoints; end

        waveforms = spikes.waveforms(spikeIDs,i1:i2,iChan);
        waveforms = psr_int16_to_single(waveforms,parameters);
        
        if (d ~= 0)
            nSpikes = size(waveforms,1);
            padding = NaN(nSpikes,abs(d));
            if     (d < 0); waveforms = [padding,waveforms];
            elseif (d > 0); waveforms = [waveforms,padding];
            end
        end

        waveformsNew(spikeIDs,:,iChan) = waveforms; 
        
    end
    
end

for iChan = 1:nChans
   % Add Gaussian noise
    waveforms = waveformsNew(:,:,iChan);
    I = find(isnan(waveforms));
    waveforms(I) = normrnd(0,spikes.info.bgn(iChan),size(I));
    waveformsNew(:,:,iChan) = waveforms; 
    
end

spikes.waveforms = int16(precision * waveformsNew);

end