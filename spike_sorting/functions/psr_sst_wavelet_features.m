function spikes = psr_sst_wavelet_features(spikes,parameters)

% From: WaveClus
% https://github.com/csn-le/wave_clus
% Unsupervised spike detection and sorting with wavelets and superparamagnetic clustering
% R. Quian Quiroga, Z. Nadasdy and Y. Ben-Shaul
% Neural Computation 16, 1661-1687; 2004.

% Input: 
%  waves - Nspikes X Npoints 

% Calculates the spike features

dims   = parameters.cluster.wavelet.dims;
scales = parameters.cluster.wavelet.scales;

waveforms = psr_int16_to_single(spikes.waveforms(:,:),parameters); 
nSpikes   = size(waveforms,1);
nSamples  = size(waveforms,2);

cc = zeros(nSpikes,nSamples);
for i = 1 : nSpikes  % Wavelet decomposition
    [c,~] = wavedec(waveforms(i,:), scales, 'haar');
    cc(i,1:nSamples) = c(1:nSamples);
end

for i = 1 : nSamples  % KS test for coefficient selection

    thr_dist     =  std(cc(:,i)) * 3;
    thr_dist_min = mean(cc(:,i)) - thr_dist;
    thr_dist_max = mean(cc(:,i)) + thr_dist;

    aux = cc(cc(:,i) > thr_dist_min & cc(:,i) < thr_dist_max, i);

    if length(aux) > 10
        [ksstat] = test_ks(aux);
        sd(i)    = ksstat; %#ok
    else
        sd(i)    = 0; %#ok
    end
end

[~, ind]      = sort(sd);
coeff(1:dims) = ind(nSamples : -1 : nSamples - dims + 1);
        
features = zeros(nSpikes,dims);
 
for i = 1 : nSpikes
    for j = 1 : dims; features(i,j) = cc(i,coeff(j)); end
end

spikes.features = single(features');

end

function [KSmax] = test_ks(x)

% Calculates the CDF (expcdf)
% [y_expcdf,x_expcdf]=cdfcalc(x);

x = x(~isnan(x));
n = length(x);
x = sort(x(:));

yCDF = (1:n)' / n; % Get cumulative sums

% Remove duplicates; only need final one with total count

notdup      = ([diff(x(:)); 1] > 0);
x_expcdf    = x(notdup);
y_expcdf    = [0; yCDF(notdup)];

% The theoretical CDF (theocdf) is assumed to be normal  
% with unknown mean and sigma

zScores = (x_expcdf - mean(x))./std(x);

mu      = 0; 
sigma   = 1; 
theocdf = 0.5 * erfc(-(zScores-mu)./(sqrt(2)*sigma));

% Compute the Maximum distance: max|S(x) - theocdf(x)|.

delta1    =  y_expcdf(1:end-1) - theocdf;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  y_expcdf(2:end)   - theocdf;   % Vertical difference at jumps approaching from the RIGHT.
deltacdf  =  abs([delta1 ; delta2]);

KSmax =  max(deltacdf);

end