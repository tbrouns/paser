function [artifacts,artifactTimes,correlations,nsamples] = ss_artifact_correlation(spikes,data,spikeTimes,criterion)

Fs = spikes.params.Fs;
artifact_length = floor(0.5 * (spikes.params.artifacts.length / 1000) * Fs); % in samples

artifactSamples = ss_artifact_filter(spikes,spikeTimes,criterion);

artifactSamples = round(Fs * artifactSamples);
artifactSamples(artifactSamples <= artifact_length) = artifact_length + 1;
artifactSamples(artifactSamples > size(data,2) - artifact_length) = size(data,2) - artifact_length;
nspikes         = length(artifactSamples);
artifactSamples = bsxfun(@plus,artifactSamples,-artifact_length:artifact_length)';
artifactSamples = artifactSamples(:);
waveforms       = data(:,artifactSamples)';
nsamples        = length(-artifact_length:artifact_length);
n               = 1;
correlations    = ones(nspikes,1);

% Calculate correlation

if (size(waveforms,2) > 1)
    for iSpike = 1:nspikes
        R = triu(corr(waveforms(n:n+nsamples-1,:)),1); % where the columns of A represent random variables and the rows represent observations.
        R = R(triu(true(size(R)),1));
        correlations(iSpike) = mean(R);
        n = n + nsamples;
    end
end

artifactSamples = reshape(artifactSamples,nsamples,nspikes);
id              = correlations > spikes.params.artifacts.corr;
artifacts       = artifactSamples(:,id);
correlations    = correlations(id);
artifactTimes   = artifacts(artifact_length + 1,:) / Fs;
artifacts       = artifacts(:);


end