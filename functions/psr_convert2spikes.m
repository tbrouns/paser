function spikes = psr_convert2spikes(spikes,data,spikepoints,assigns,parameters)

% PSR_CONVERT2SPIKES - Convert spike times and cluster assigns to PASER data format.
%
% Syntax:  spikes = psr_convert2spikes(rez,data,parameters)
%
% Inputs:
%    spikes     - Can be empty 
%    data       - Filtered time series of extracellular recording [Nchannels x Npoints]
%    spiketimes - Sample points of each spike [Nspikes x 1] 
%    assigns    - Cluster index of each spike [Nspikes x 1]
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

%% Save results

spikes.assigns    = int16(assigns');
spikes.spiketimes = single((spikepoints - 1)' / Fs); % in secs
spikes.waveforms  = psr_sst_get_waveforms(spikepoints,data,win);

%------------- END OF CODE --------------