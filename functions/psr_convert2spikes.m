function spikes = psr_convert2spikes(data,spiketimes,assigns,parameters)

% PSR_CONVERT2SPIKES - Convert spike times and cluster assigns to PASER data format.
%
% Syntax:  spikes = psr_convert2spikes(rez,data,parameters)
%
% Inputs:
%    data - Filtered time series of extracellular recording [number of
%    probe channels x number of data points]
%    spiketimes - Time in seconds of each spike [number of spikes x 1] 
%    assigns - Cluster index of each spike [number of spikes x 1]
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
sLength = size(data,2);
    
% Get raw data around spiketimes

times = bsxfun(@plus,spiketimes,win);
times(times < 1) = 1;
times(times > sLength) = sLength;
times = times';
times = times(:);
waves = data(:,times);
waves = permute(waves,[3 2 1]);
waves = reshape(waves,length(win),[],parameters.general.nelectrodes);
waves = permute(waves,[2 1 3]);

%% Save results

spikes = [];
spikes.assigns     = int16(assigns');
spikes.spiketimes  = single((spiketimes - 1)' / Fs); % in secs
spikes.waveforms   = waves;

%------------- END OF CODE --------------