function [artifacts,artifactTimes,correlations,nsamples] = psr_sst_artifact_correlation(parameters,data,spikeTimes)

Fs = parameters.Fs;
artifact_length = floor(0.5*Fs*(parameters.spikes.window_size / 1000)); % in samples

[artifactSamples,artifactIDs] = psr_sst_align(parameters,spikeTimes); % align adjacent spikes

artifactSamples = round(Fs * artifactSamples);
artifactSamples(artifactSamples <= artifact_length) = artifact_length + 1;
artifactSamples(artifactSamples > size(data,2) - artifact_length) = size(data,2) - artifact_length;
nspikes      = length(artifactSamples);
nsamples     = length(-artifact_length:artifact_length);
correlations = ones(nspikes,1);

% Calculate correlation

for iSpike = 1:nspikes

    samples = artifactSamples(iSpike);
    samples = bsxfun(@plus,samples,-artifact_length:artifact_length)';
    samples = samples(:);
    waveforms = data(:,samples)';

    R = triu(corr(waveforms),1); % where the columns of A represent random variables and the rows represent observations.
    R = R(triu(true(size(R)),1));
    correlations(iSpike) = mean(R);

end

id              = correlations > parameters.spikes.artifacts_corr;
artifacts       = artifactSamples(id);
correlations    = correlations(artifactIDs);
artifactTimes   = artifacts / Fs;

end