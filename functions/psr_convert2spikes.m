function spikes = psr_convert2spikes(spikes,data,spikepoints,assigns,parameters)

% PSR_CONVERT2SPIKES - Convert to PASER spike data format 
% This function takes spike sample points and cluster assigns as inputs and
% generates a structure that complies with the PASER data format
%  
% Syntax:  spikes = psr_convert2spikes(spikes,data,spikepoints,assigns,parameters)
%
% Inputs:
%    spikes      - Structure to add data to. Can also be empty.
% 
%    data        - Matrix of band-pass filtered voltage time series, with
%                  shape:
%                  [Number of channels x Number of data points]
% 
%    spikepoints - Vector of sample points of each spike, with shape:
%                  [Number of spikes x 1] 
% 
%    assigns     - Vector of cluster indices of all spikes, with shape:
%                  [Number of spikes x 1]
% 
%    parameters  - See README
%
% Outputs:
%    spikes - See README
% 
% See also: PSR_SST_SORTING_KST,  PSR_SST_GET_WAVEFORMS

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