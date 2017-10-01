function [spikeTimesAligned,spikeIDs] = ept_sst_align(parameters,spikeTimes)

% Filter artifacts that do not occur on enough channels

spikeTimes        =   sort(spikeTimes);
nspikes           = length(spikeTimes);
spikeDur          = (parameters.spikes.artifacts_corr * parameters.spikes.window_size) / 1000; % sec
spikeTimesAligned = -1 * ones(nspikes,1);
spikeIDs          = zeros(nspikes,1);
iSpike            = 1;
kSpike            = 1;
T                 = spikeTimes(1);

while (iSpike <= nspikes)
    spikeTimesTemp = [];
    jSpike = 1;
    while (iSpike <= nspikes) % end condition
        t = spikeTimes(iSpike);
        dt = t - T; T = t;
        if (dt > spikeDur); break; end
        spikeTimesTemp(jSpike) = t; %#ok
        spikeIDs(iSpike) = kSpike;
        iSpike = iSpike + 1;        
        jSpike = jSpike + 1;
    end
    spikeTimesAligned(kSpike) = mean(spikeTimesTemp);
    kSpike = kSpike + 1;
end

spikeTimesAligned(spikeTimesAligned < 0) = [];

end