function [x,n] = psr_sst_amp_hist(spikes,clustID,parameters,splitAmp)

if (splitAmp); amplitudes = psr_sst_amp_split(spikes,clustID,parameters);
else,          amplitudes = psr_sst_amp      (spikes,clustID,parameters);
end

th(1,1,:) = mean(spikes.info.thresh);
minBin = 10 * (10^-parameters.general.precision) / abs(mean(th(:)));
bin = range(amplitudes) * (parameters.cluster.amplitude_nbin / length(amplitudes));
if (bin < minBin); bin = minBin; end
X = min(amplitudes):bin:max(amplitudes);
n = histcounts(amplitudes,X);
x = X(1:end-1) + (X(2) - X(1))/2;
if all(spikes.info.thresh < 0); x = -x; end

end
