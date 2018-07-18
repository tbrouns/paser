function spikes = psr_freq2spikes(spikes,freq)

if (isfield(freq,'artifacts'))
    nTrials = size(freq,2);
    for iTrial = 1:nTrials
        x = freq(iTrial).artifacts;
        fields = fieldnames(x);
        nFields = size(fields,1);
        artifacts = [];
        for iField = 1:nFields
            artifacts.(fields{iField}) = (freq(iTrial).artifacts.(fields{iField}) - 1) / freq(iTrial).fsample; 
        end
        spikes.info.artifacts.lfp{iTrial,1} = artifacts;
    end
end
end