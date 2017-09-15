function nlength = ss_artifact_removal(files,data,MFAtimes_1,MFAtimes_2)

% TO DO: Change parameters to default spike parameters

% Artifact removal

spikeTimes = [];
ntets      = length(files);
Fs         = zeros(ntets,1);

for iTetrode = 1:ntets
    load(files{iTetrode});
    nlength = spikes.nlength; %#ok
    nmax = max(nlength);
    if (isfield(spikes,'spiketimes'))
        for i = 1:nmax
            spikeTimes = [spikeTimes, spikes.spiketimes(find(nlength))]; %#ok
            nlength = nlength - 1;
        end
    end
    Fs(iTetrode) = spikes.params.Fs;
end

Fs = mode(Fs);

% Only consider artifacts that occur on at least 2 channels simultaneously

[artifacts,~,~,nsamples] = ss_artifact_correlation(spikes,data,spikeTimes,2);

% Add magnetic field artifacts

MFAtimes = [MFAtimes_1, MFAtimes_2];
MFAtimes = MFAtimes';

if (~isempty(MFAtimes))
    artifact_length = 0.5 * (nsamples - 1);
    MFAs = round(Fs * MFAtimes);
    MFAs(MFAs <= artifact_length) = artifact_length + 1;
    MFAs(MFAs > size(data,2) - artifact_length) = size(data,2) - artifact_length;
    MFAs = bsxfun(@plus,MFAs,-artifact_length:artifact_length)';
    MFAs = MFAs(:);
    artifacts = [artifacts; MFAs];
end

% Process

nlength = zeros(ntets,1);
for iTetrode = 1:ntets
    
    samples = zeros(1,size(data,2));
    load(files{iTetrode});
        
    % Save correlations, artifact times and waveforms
    
    spikes.artifacts_1 = single(MFAtimes_1);
    spikes.artifacts_2 = single(MFAtimes_2);
    
    %     spikes.artifact_correlations = single(correlations');
    
    % Remove artifact spikes from spike vector
    
    spiketimes = spikes.spiketimes; % in seconds
    spiketimes = round(spiketimes * Fs);   % in sample number
    samples(spiketimes)  = 2;
    samples(artifacts)   = samples(artifacts) - 1;
    samples(samples < 1) = [];
    id = find(samples == 2);
    spikes = ss_spike_removal(spikes,id,0);
    
    % Save
    save(files{iTetrode},'spikes','freq','parameters');
    nlength(iTetrode) = length(id);
end

end