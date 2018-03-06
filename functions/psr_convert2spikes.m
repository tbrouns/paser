function spikes = psr_convert2spikes(spikes,data,spiketimes,assigns,parameters)

% PSR_CONVERT2SPIKES - Convert spike times and cluster assigns to PASER data format.
%
% Syntax:  spikes = psr_convert2spikes(rez,data,parameters)
%
% Inputs:
%    spikes - Can be empty 
%    data - Filtered time series of extracellular recording [Nchannels x Npoints]
%    spiketimes - Data points of each spike [Nspikes x 1] 
%    assigns - Cluster index of each spike [Nspikes x 1]
%    parameters - See README
%
% Outputs:
%    spikes - See README

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

Fs = parameters.Fs;
window_samples = round(Fs * parameters.spikes.window_size / 1000);
samples_hwidth = round(0.5 * window_samples);
win = -samples_hwidth:samples_hwidth;

% Get raw data around spiketimes

waveforms = psr_sst_get_waveforms(spiketimes,data,win);

% sLength = size(data,2);
% timeArray = bsxfun(@plus,spiketimes,win);
% timeArray(timeArray < 1) = 1;
% timeArray(timeArray > sLength) = sLength;
% timeArray = timeArray';
% timeArray = timeArray(:);
% waveforms = data(:,timeArray);
% waveforms = permute(waveforms,[3 2 1]);
% waveforms = reshape(waveforms,length(win),[],parameters.general.nelectrodes);
% waveforms = permute(waveforms,[2 1 3]);

%% Save results

spikes.assigns    = int16(assigns');
spikes.spiketimes = single((spiketimes - 1)' / Fs); % in secs
spikes.waveforms  = waveforms;

%------------- END OF CODE --------------