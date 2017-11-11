function spikes = psr_sst_detect(signal, spikes, parameters)

% psr_sst_detect - Detection of spikes in multi-channel data
%
% Adapted from 'ss_detect' originally created for UltraMegaSort2000 by Hill
% DN, Mehta SB & Kleinfeld D (2010). Major edits made by Terence Brouns
% (2017).
% 
% Description:
%     SPIKES = PSR_SST_DETECT(data, spikes) takes a matrix of data and
%     returns a spikes object containing the waveforms and spike times
%     embedded in that signal.  This is determined by observing threshold
%     crossings on all channels. By default, the threshold is determined by
%     maximum absolute deviation (MAD).  All threshold crossings on any
%     channel are events as long as they are not preceeded too closely be
%     another event.  See SPIKES.PARAMS.SHADOW below.
%
%     NOTE: it is assumed that the data has been previously filtered and
%     that there is no voltage offset to the signal so that the mean signal
%      is 0.
%
%     Input data must be in a matrix format [samples X channels].
%     If a data set is too large to fit in memory, ss_detect can be called
%     multiple times with different data sets.  The new spikes will be
%     appended to the data for the previously detected spikes.
%
%     This function saves a window of data samples from each channel for
%     each detected event.  This window is determined by parameters in
%     spikes.params and described below.  In brief, the user specifies
%     (1) how long a window to extract, (2) where in the window the
%     threshold crossing shuld appear, and (3) a minimum separation time
%     between detected events.
%
%     The function places the following fields in the spikes object:
%       SPIKES.WAVEFORMS   : [events X SAMPLES X CHANNELS] matrix of event waveforms
%       SPIKES.SPIKETIMES  : array containing time within trial of each event
%       SPIKES.INFO.DETECT : structure containing information on detect session - see manual for more details
%
%     TODO: add separate parameter structure

append = isfield(spikes, 'waveforms');

if ~append
    spikes.waveforms  = [];
    spikes.spiketimes = [];
    spikes.nspikes    = [];
    spikes.info.dur = 0;
end

% set some constants
Fs          = spikes.Fs;
method      = parameters.spikes.method;
sLength     = size(signal,1); % total number of samples in signal
nChan       = size(signal,2);
sWindow     = round(Fs * parameters.spikes.window_size / 1000);
sWindowHalf = round(0.5 * sWindow);

% determine threshold

stdev = std(signal);
if     isequal(method, 'auto');   thresh = -parameters.spikes.thresh * stdev;
elseif isequal(method, 'manual'); thresh =  parameters.spikes.thresh;
elseif isequal(method, 'mad');    thresh = -parameters.spikes.thresh * psr_mad(signal);
else,  error('Unknown spike detection method.')
end

spikes.info.stds   = stdev;
spikes.info.thresh = thresh;
    
% Get crossings on all channels for this trial

crossings = zeros(nChan,sLength);
peaks     = zeros(nChan,sLength);

for k = 1:nChan
    signalChan = signal(:,k);
    del = signalChan > thresh(k);
    I   = find(~del); % keep raw locations of thresholded data
    signalChan(del) = []; % ignore sub-threshold data
    if (length(signalChan) > 3) % condition required for 'findpeaks'
        [pks,loc] = findpeaks(double(sign(thresh(k)) * signalChan));
        loc = I(loc); % find locations of peaks in raw data
        crossings(k,loc) = 1;
        peaks(k,loc)     = pks;
    end
end

divisor   = abs(repmat(spikes.info.thresh', 1, size(peaks,2)));
peaks     = peaks ./ divisor;
peaksMax  = max(peaks);
crossings = find(sum(crossings));
peaksMax  = peaksMax(crossings);

% Deal with adjacent peaks, within spike window. 
% Keep largest amplitude peak

indices = diff(crossings) <= sWindowHalf;
indices = [indices,0];
nspikes = length(indices);
indices = [indices,0]; % needed for stop criterion
peakCrossings = [];
iSpike = 1;
while iSpike <= nspikes
    pks = [];
    itr = 0;
    while(itr == 0 || indices(iSpike+itr-1))
        pks = [pks;peaksMax(iSpike+itr)]; %#ok
        itr = itr + 1;
    end
    [~,I] = max(pks);
    peakCrossings = [peakCrossings;crossings(iSpike+I-1)]; %#ok
    iSpike = iSpike + itr;
end

crossings = peakCrossings';
crossings(crossings <= sWindowHalf) = [];
crossings(crossings > sLength - sWindowHalf) = [];

% update spiketimes, trials, and waveforms
spikes.spiketimes = [spikes.spiketimes, (crossings - 1) / Fs];
indices = bsxfun(@plus,crossings,-sWindowHalf:sWindowHalf); 
waves   = signal(indices, :);

spikes.waveforms = [spikes.waveforms; waves];

% Duration of trial
    
clear data

% Save everything
spikes.waveforms  = single(spikes.waveforms);
spikes.spiketimes = single(spikes.spiketimes);

% report detection rate
detect_rate = length(spikes.spiketimes) / sum(spikes.info.dur);
disp(['Detected on average ' num2str(detect_rate) ' events per second of data ']);