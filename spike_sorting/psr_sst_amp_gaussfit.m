function [p,mu,stdev] = psr_sst_amp_gaussfit(spikes,clustID,parameters)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_undetected - wrapper for undetected.m
%
% Usage:
%     [p,mu,stdev,criteria] = ss_undetected(spikes,use)
%
% Description:
%    Estimates how fraction of spikes associated with a given cluster that
% may have gone undetected, by fitting a Gaussian with a missing tail to a
% histogram of the detection criterion applied to each spike event.
%
% See undetected.m for more information.
%
% Input:
%   spikes - spikes structure
%   clustID - a cluster ID or an array describing which spikes to
%            use as the first cluster
%
% Output:
%  p            - estimate of probability that a spike is missing because it didn't reach threshhold
% 
% 
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% undetected - estimate fraction of events that did not reach threshold
%
% Usage:
%     [p,mu,stdev,n,x] = undetected(waveforms,threshes,criteria_func)
%
% Description:
%    Estimates fraction of events that did not reach threshold by applying
% the detection metric to each waveform and then fitting it with a Gaussian
% that has a missing tail.
%
% The distribution of detection metric values is turned into a histogram.
% A Gaussian is fit to the historgram to minimize the absolute error
% between the Gaussian and the histogram for values above threshold.
% The integral of this Gaussian that is below threshold is the estimate of
% the fraction of missing events.
%
% Note that values are normalized so that the threshold is +/- 1.  The function
% attempts to preserve the sign of the original threshold, unless thresholds
% on different channels had different signs. In the case of multiple channels,
% each channel is normalized so that the threshold has a magnitude of 1.  Then,
% for each event, only the channel with the most extreme value of the detection
% metric is used.

% create the histogram values

mu    = [];
stdev = [];

if (spikes.info.detected) % If spike threshold was used
    
    [x,n,amplitudes] = psr_sst_amp_hist(spikes,clustID,parameters,true);
    
    % fit the histogram with a cutoff gaussian
    m = mode_guesser(amplitudes, 0.05);    % use mode instead of mean, since tail might be cut off
    [stdev,mu] = stdev_guesser(amplitudes, n, x, m); % fit the standard deviation as well

    % Now make an estimate of how many spikes are missing, given the Gaussian and the cutoff
    p = normcdf(1,mu,stdev);
    
    % attempt to keep values negative if all threshold values were negative
    if all(spikes.info.thresh < 0); mu = -mu; end

else
    spikeIDs  = ismember(spikes.assigns,clustID);
    waveforms = spikes.waveforms(spikeIDs,:,:);
    waveforms = psr_int16_to_single(waveforms,parameters);
    waveforms = psr_sst_norm_waveforms(waveforms,mean(spikes.info.thresh));
    waveforms = max(waveforms,[],2);
    waveforms = max(waveforms,[],3);
    nspikes = size(waveforms,1);
    p = sum(waveforms < 1) / nspikes;
end

end

% fit the standard deviation to the histogram by looking for an accurate
% match over a range of possible values
function [stdev,m] = stdev_guesser(thresh_val, n, x, m)

% initial guess is juts the RMS of just the values below the mean
init = sqrt(mean((m - thresh_val(thresh_val >= m)).^2));

% try 20 values, within a factor of 2 of the initial guess
num        = 20;
sd_guesses = linspace(init/3, init*3, num);
Nsd = length(sd_guesses);
errors = zeros(1,Nsd);
for iStd = 1:Nsd
    b = normpdf(x,m,sd_guesses(iStd));
    b = b * max(n) / max(b);
    errors(iStd) = sum(abs(b(:)-n(:)));
end

% which one has the least error?
[~,pos] = min(errors(:));
kpos    = ceil(pos / num);
stdev   = sd_guesses(kpos);


end