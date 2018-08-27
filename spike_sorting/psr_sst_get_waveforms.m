function waveforms = psr_sst_get_waveforms(spikePoints,data,win)

% data - Filtered time series of extracellular recording [Nchannels x Npoints]
% spiketimes - Sample number of each spike [Nspikes x 1] 

nChans  = size(data,1);
sLength = size(data,2);
spikePoints = double(spikePoints); % Avoids erroneous results
spikePoints = spikePoints(spikePoints <= sLength);
timeArray = bsxfun(@plus,spikePoints,win);
timeArray(timeArray < 1) = 1;
timeArray(timeArray > sLength) = sLength;
timeArray = timeArray';
timeArray = timeArray(:);
waveforms = data(:,timeArray);
waveforms = permute(waveforms,[3 2 1]);
waveforms = reshape(waveforms,size(win,2),[],nChans);
waveforms = permute(waveforms,[2 1 3]);

end