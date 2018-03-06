function waveforms = psr_sst_get_waveforms(spiketimes,data,win)

% data - Filtered time series of extracellular recording [Nchannels x Npoints]
% spiketimes - Time in seconds of each spike [Nspikes x 1] 

nChans  = size(data,1);
sLength = size(data,2);
timeArray = bsxfun(@plus,spiketimes,win);
timeArray(timeArray < 1) = 1;
timeArray(timeArray > sLength) = sLength;
timeArray = timeArray';
timeArray = timeArray(:);
waveforms = data(:,timeArray);
waveforms = permute(waveforms,[3 2 1]);
waveforms = reshape(waveforms,length(win),[],nChans);
waveforms = permute(waveforms,[2 1 3]);

end