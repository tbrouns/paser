function data = psr_bp_filter(data,cfg)

% PSR_BP_FILTER - Zero-phase digital band-pass filtering
%
% Syntax:  data = psr_bp_filter(data,cfg)
%
% Inputs:
%    data - See input signal for FILTFILT function 
%    cfg  - Structure containing the following fields:
%           "order" : Order of Butterworth filter
%           "lower" : Lower  cut-off frequency
%           "upper" : Higher cut-off frequency
%           "Fs"    : Sampling frequency of input signal
%
% Outputs:
%    data - Band-pass filtered signal
% 
% Dependencies: Signal Processing Toolbox
%
% See also: FILTFILT,  BUTTER

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

[B,A] = butter(cfg.order,[cfg.lower cfg.upper]/(cfg.Fs/2),'bandpass');
data  = filtfilt(B,A,data); % Zero-phase digital filtering

end