function spikesMain = psr_sst_spike_append(spikesMain,spikes)

spiketimes = spikes.spiketimes;
tstart     = sum(spikesMain.info.dur);
if (isempty(tstart)); tstart = 0; end
spiketimes = spiketimes + tstart;

spikesMain.waveforms  = [spikesMain.waveforms; spikes.waveforms];
spikesMain.spiketimes = [spikesMain.spiketimes,spiketimes];

if (isfield(spikesMain,'assigns') && isfield(spikes,'assigns'))
    spikesMain.assigns = [spikesMain.assigns, spikes.assigns];
end

if (isfield(spikesMain,'data') && isfield(spikes,'data'))
    spikesMain.data = [spikesMain.data,spikes.data];
end

if (isfield(spikesMain,'info') && isfield(spikes,'info'))
    
    if (isfield(spikesMain.info,'detect') && isfield(spikes.info,'detect'))
        
        if (isfield(spikesMain.info,'stds') && isfield(spikes.info,'stds'))
            spikesMain.info.stds   = [spikesMain.info.stds   ; spikes.info.stds];
        end
        
        if (isfield(spikesMain.info,'thresh') && isfield(spikes.info,'thresh'))
            spikesMain.info.thresh = [spikesMain.info.thresh ; spikes.info.thresh];
        end
        
        if (isfield(spikesMain.info,'dur') && isfield(spikes.info,'dur'))
            spikesMain.info.dur    = [spikesMain.info.dur    ; spikes.info.dur];
        end
        
    end
    
end

end