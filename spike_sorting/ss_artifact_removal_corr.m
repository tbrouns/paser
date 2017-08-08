function nspikes = ss_artifact_removal_corr(files, data)

% Artifact removal
    
spikeTimes = [];
ntets      = length(files);

for iTetrode = 1:ntets
    load(files{iTetrode});
    spikeTimes = [spikeTimes; spikes.artifacts]; %#ok
end

Fs = spikes.params.Fs;

artifact_length = floor(0.5 * (spikes.params.artifact_length / 1000) * spikes.params.Fs); % in samples

% Find tentative artifacts

spikeTimes    = sort(spikeTimes);
nspikes       = length(spikeTimes);
criterion     = round(ntets * spikes.params.artifact_fract);
if (criterion < 2); criterion = 2; end % to avoid any errors
artifactDur   = round(Fs * (spikes.params.artifact_offset / 1000));
artifactSamples = -1 * ones(nspikes,1);
iSpike        = 1;
kSpike        = 1;

while iSpike < nspikes
    jSpike = iSpike + 1;
    t      = spikeTimes(iSpike);
    dt     = 0;
    while(dt <= artifactDur && jSpike <= nspikes)
        dt = spikeTimes(jSpike) - t;
        jSpike = jSpike + 1;
    end
    if (jSpike - iSpike - 1 > criterion)
        artifactSamples(kSpike) = 0.5 * (spikeTimes(jSpike - 2) + t);
        kSpike = kSpike + 1;
    end
    iSpike = jSpike - 1;
end

artifactSamples(artifactSamples < 0) = [];
artifactSamples = round(artifactSamples);
artifactSamples(artifactSamples <= artifact_length) = artifact_length + 1;
artifactSamples(artifactSamples > size(data,2) - artifact_length) = size(data,2) - artifact_length;
nspikes         = length(artifactSamples);
artifactSamples = bsxfun(@plus,artifactSamples,-artifact_length:artifact_length)';
artifactSamples = artifactSamples(:);
waveforms       = data(:,artifactSamples)';
nsamples        = length(-artifact_length:artifact_length);
n               = 1;
correlations    = zeros(nspikes,1);

for iSpike = 1:nspikes
    R = triu(corr(waveforms(n:n+nsamples-1,:)),1); % where the columns of A represent random variables and the rows represent observations.
    R = R(triu(true(size(R)),1));
    correlations(iSpike) = mean(R);
    n = n + nsamples;
end

artifactSamples = reshape(artifactSamples,nsamples,nspikes);
id              = correlations > spikes.params.artifact_corr;
artifacts       = artifactSamples(:,id);
correlations    = correlations(id);
artifactTimes   = artifacts(artifact_length + 1,:) / Fs;
artifacts       = artifacts(:);

nspikes = zeros(ntets,1);
for iTetrode = 1:ntets
    samples = zeros(1,size(data,2));
    load(files{iTetrode});
    spikes.artifacts    = single(artifactTimes);
    spikes.correlations = single(correlations');
    spiketimes = spikes.spiketimes; % in seconds
    spiketimes = round(spiketimes * Fs);   % in sample number
    samples(spiketimes)  = 2;
    samples(artifacts)   = samples(artifacts) - 1;
    samples(samples < 1) = []; 
    id = find(samples == 2);
    spikes.spiketimes      = spikes.spiketimes(id);
    spikes.trials          = spikes.trials(id);
    spikes.waveforms       = spikes.waveforms(id,:,:);
    spikes.unwrapped_times = spikes.unwrapped_times(id);
    save(files{iTetrode},'spikes');
    nspikes(iTetrode) = length(id);
end

end