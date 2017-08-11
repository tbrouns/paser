function artifactSamples = ss_artifact_filter(spikes,spikeTimes,criterion)

spikeTimes    =   sort(spikeTimes);
nspikes       = length(spikeTimes);
artifactDur   = spikes.params.artifact_offset / 1000;
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

end