function signal = psr_stability_reversal(spikes,signal,parameters)

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);

Fs = parameters.Fs;
window_samples = round(Fs * parameters.spikes.window_size / 1000);
samples_hwidth = round(0.5 * window_samples);
win = -samples_hwidth:samples_hwidth;

for iClust = 1:nClust
    clusterID = clusterIDs(iClust);
    which     = (spikes.assigns == clusterID);
    t = spikes.spiketimes(which);
    nspikes = length(t);
    t = (t / Fs) + 1;
    t = bsxfun(@plus,t,win);
    wmean = signal(t);
    wmean = mean(wmean,1);
    wmean = repmat(wmean,nspikes,1);
    waveSignal = zeros(size(signal));
    waveSignal(t) = wmean;
    signal = 2 * waveSignal(t) - signal;
end

end