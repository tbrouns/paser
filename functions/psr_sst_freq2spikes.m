function spikes = psr_sst_freq2spikes(spikes,freq)

if (isfield(freq,'artifacts'))
    nTrials = size(freq,2);
    for iTrial = 1:nTrials
        artifacts = (freq(iTrial).artifacts - 1) / freq(iTrial).fsample; 
        spikes.info.artifacts{iTrial,1} = artifacts;
    end
end
end