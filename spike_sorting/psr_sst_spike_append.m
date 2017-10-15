function spikesMain = psr_sst_spike_append(spikesMain,spikes)

spiketimes = spikes.spiketimes;
tstart     = sum(spikesMain.info.detect.dur);
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
        
        if (isfield(spikesMain.info.detect,'stds') && isfield(spikes.info.detect,'stds'))
            spikesMain.info.detect.stds   = [spikesMain.info.detect.stds   ; spikes.info.detect.stds];
        end
        
        if (isfield(spikesMain.info.detect,'thresh') && isfield(spikes.info.detect,'thresh'))
            spikesMain.info.detect.thresh = [spikesMain.info.detect.thresh ; spikes.info.detect.thresh];
        end
        
        if (isfield(spikesMain.info.detect,'dur') && isfield(spikes.info.detect,'dur'))
            spikesMain.info.detect.dur    = [spikesMain.info.detect.dur    ; spikes.info.detect.dur];
        end
        
    end
    
end

end