function nspikes = ss_artifact_removal(files,data)

% TO DO: Change parameters to default spike parameters

% Artifact removal

MFAtimes   = [];
MFAs       = [];
spikeTimes = [];
ntets      = length(files);
nchan      = 4;

for iTetrode = 1:ntets
    load(files{iTetrode});
    nspikes = spikes.nspikes;
    nmax = max(nspikes);
    for i = 1:nmax
        spikeTimes = [spikeTimes, spikes.spiketimes(find(nspikes))]; %#ok
        nspikes = nspikes - 1;
    end
    
    if isfield(spikes,'artifacts')
       MFAtimes = [MFAtimes, spikes.artifacts]; %#ok 
    end
end

Fs = spikes.params.Fs;

% Filter artifacts based on correlation across all channels

criterion     = round(4 * ntets * spikes.params.artifact_fract);
if (criterion < 2); criterion = 2; end % to avoid any errors

[artifacts,~,correlations,nsamples] = ss_artifact_correlation(spikes,data,spikeTimes,criterion);

% Add magnetic field artifacts

if (~isempty(MFAtimes))
    artifact_length = 0.5 * (nsamples - 1);
    MFAtimes = ss_artifact_filter(spikes,MFAtimes,0);
    MFAs     = MFAtimes;
    MFAtimes = round(Fs * MFAtimes);
    MFAtimes(MFAtimes <= artifact_length) = artifact_length + 1;
    MFAtimes(MFAtimes > size(data,2) - artifact_length) = size(data,2) - artifact_length;
    MFAtimes = bsxfun(@plus,MFAtimes,-artifact_length:artifact_length)';
    MFAtimes = MFAtimes(:);
    artifacts = [artifacts; MFAtimes];
end

% Process

nspikes = zeros(ntets,1);
for iTetrode = 1:ntets
    
    samples = zeros(1,size(data,2));
    load(files{iTetrode});
    
    % Save correlations, artifact times and waveforms
    if (~isempty(MFAs)); spikes.artifacts = single(MFAs); end
    waveforms = data(nchan*(iTetrode-1)+1:nchan*iTetrode,artifacts);
    waveforms = permute(waveforms,[2 3 1]);
    waveforms = reshape(waveforms,nsamples,[],nchan);
    spikes.artifact_waveforms = single(permute(waveforms,[2 1 3]));
    spikes.artifact_correlations = single(correlations');
    
    % Remove artifact spikes from spike vector
        
    spiketimes = spikes.spiketimes; % in seconds
    spiketimes = round(spiketimes * Fs);   % in sample number
    samples(spiketimes)  = 2;
    samples(artifacts)   = samples(artifacts) - 1;
    samples(samples < 1) = []; 
    id = find(samples == 2);
    spikes = ss_spike_removal(spikes,id,0);
    
    % Save
    save(files{iTetrode},'spikes');
    nspikes(iTetrode) = length(id);
end

end