function spikes = psr_sst_background_noise(spikes,data,parameters)

% Calculate background noise of each channel and duration of each trial
%
% 'env': Background noise level is calculated using method presented in [1]
%
% References:
% [1] Dolan, Kevin, et al. "Automatic noise-level detection for extra-cellular
% micro-electrode recordings." Medical & biological engineering & computing
% 47.7 (2009): 791-800.

nChans = size(data,1);

for iChan = 1:nChans
    
    dataChan = data(iChan,:);    
    dataChan = psr_int16_to_single(dataChan,parameters);
    spikes.info.std(1,iChan) =     std(dataChan,[],2)'; % Standard deviation
    spikes.info.mad(1,iChan) = psr_mad(dataChan     )'; % Median absolute deviation
    spikes.info.rms(1,iChan) =     rms(dataChan,   2)'; % Root-mean-square
    
    % Find mode of signal envelope: magnitude of analytical signal
    
    [yUpper,~] = envelope(dataChan);

    binsize = 10^-parameters.general.precision;
    binmax  = mean(yUpper) + 5 * std(yUpper);
    edges   = 0:binsize:binmax;
    [count,edges] = histcounts(yUpper,edges);
    edges = edges(1:end-1) + 0.5 * diff(edges);
    [~,I] = max(count);
    spikes.info.env(1,iChan) = edges(I);

end

% Set background noise
switch parameters.spikes.bgntype % Background noise type
    case 'env'; spikes.info.bgn = spikes.info.env;
    case 'mad'; spikes.info.bgn = spikes.info.mad;
    case 'std'; spikes.info.bgn = spikes.info.std;
    case 'rms'; spikes.info.bgn = spikes.info.rms;
end

end
