function spikes = ept_sst_detect(data, spikes, parameters)

% ept_sst_detect - Detection of spikes in multi-channel data
%
% Adapted from 'ss_detect' originally created for UltraMegaSort2000 by Hill
% DN, Mehta SB & Kleinfeld D (2010). Major edits made by Terence Brouns
% (2017).
% 
% Description:
%     SPIKES = EPT_SST_DETECT(data, spikes) takes a matrix of data and
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
    spikes.info.detect.dur = 0;
end

% set some constants
Fs             = spikes.params.Fs;
method         = parameters.spikes.method;
num_channels   = size(data, 2);
window_samples = round(Fs * parameters.spikes.window_size / 1000);
shadow         = round(Fs * parameters.spikes.shadow      / 1000);
samples_before = round(Fs * parameters.spikes.cross_time  / 1000);
samples_after  = round(Fs * parameters.spikes.max_jitter  / 1000) + window_samples - (1 + samples_before);
jitter_range   = samples_before - 1 + (1:round(parameters.spikes.max_jitter * Fs / 1000));

% determine threshold

spikes.info.detect.cov = get_covs(data, window_samples);

stdev = std(data);
if     isequal(method, 'auto');   thresh = -parameters.spikes.thresh * stdev;
elseif isequal(method, 'manual'); thresh =  parameters.spikes.thresh;
elseif isequal(method, 'mad');    thresh = -parameters.spikes.thresh * (median(abs(data)) / 0.6745);
else   error('Unknown spike detection method.')
end

spikes.info.detect.stds   = stdev;
spikes.info.detect.thresh = thresh;
    
% Get crossings on all channels for this trial

crossings = zeros(num_channels,size(data,1)-1);
peaks     = zeros(num_channels,size(data,1)-1);

for k = 1:num_channels
    data_channel = data(:,k);
    del = data_channel > thresh(k);
    I   = find(~del); % keep raw locations of thresholded data
    data_channel(del) = []; % ignore sub-threshold data
    [pks,loc] = findpeaks(double(sign(thresh(k)) * data_channel));
    loc = I(loc); % find locations of peaks in raw data
    crossings(k,loc) = 1;
    peaks(k,loc)     = pks;
end

divisor   = abs(repmat(spikes.info.detect.thresh', 1, size(peaks,2)));
peaks     = peaks ./ divisor;
peaksMax  = max(peaks);
crossings = find(sum(crossings));
peaksMax  = peaksMax(crossings);

% Deal with adjacent peaks, within shadow period. 
% Keep largest amplitude peak

indices = diff(crossings) <= shadow;
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
crossings(crossings <= samples_before) = [];
crossings(crossings > size(data,1) - samples_after) = [];

% update spiketimes, trials, and waveforms
spikes.spiketimes = [spikes.spiketimes crossings / Fs];
w = zeros([length(crossings), samples_before + 1 + samples_after, num_channels], 'single');

for k = 1:length(crossings) % maybe use bsxfun here for speed-up
    indices  = crossings(k) + (-samples_before:samples_after);
    w(k,:,:) = data(indices, :);
end
spikes.waveforms = [spikes.waveforms; w];

% get number of above threshold channels for each spike

nspikes = zeros(1,length(crossings));
for k = 1:num_channels
    peaks = max(sign(thresh(k)) * w(:,jitter_range,k),[],2);
    peaks = peaks > abs(thresh(k));
    nspikes = nspikes + peaks';
end
spikes.nspikes = [spikes.nspikes nspikes];

% Duration of trial

spikes.info.detect.dur = spikes.info.detect.dur + (size(data, 1) / Fs);
    
clear data

% Save everything
spikes.waveforms  = single(spikes.waveforms);
spikes.spiketimes = single(spikes.spiketimes);

% save some more data that will be useful later
spikes.info.detect.align_sample = samples_before + 1;

% report detection rate
detect_rate = length(spikes.spiketimes) / sum(spikes.info.detect.dur);
disp(['Detected on average ' num2str(detect_rate) ' events per second of data ']);

% get covariance matrix of background noise by randomly sampling 10000 timepoints
function c = get_covs(data, samples)

num_trials   = length(data);
num_channels = size(data{1},2);
num_samples  = zeros(1,num_trials);

for j = 1:num_trials, num_samples(j) = size(data{j},1); end

max_samples = 10000;
waves       = zeros([max_samples samples num_channels]);
tr_index    = ceil(num_trials * rand([1 max_samples]));
data_index  = ceil((num_samples(tr_index)-samples) .* rand([1 max_samples]));

for j = 1:max_samples
    waves(j,:,:) = data{tr_index(j)}(data_index(j) + (0:samples-1),:);
end

c = cov(waves(:,:));