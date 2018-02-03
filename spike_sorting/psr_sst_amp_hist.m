function [x,n,A] = psr_sst_amp_hist(spikes,clustID,parameters,SPLIT)

if (SPLIT); A = psr_sst_amp_split(spikes,clustID,parameters);
else,       A = psr_sst_amp      (spikes,clustID,parameters);
end

thresh = abs(mean(spikes.info.thresh(:)));
minBin = 10 * (10^-parameters.general.precision) / thresh;
bin = range(A) * (parameters.cluster.amplitude_nbin / length(A));
if (bin < minBin); bin = minBin; end
X = min(A):bin:max(A);
n = histcounts(A,X);
x = X(1:end-1) + (X(2) - X(1))/2;
if all(spikes.info.thresh < 0); x = -x; end

end
