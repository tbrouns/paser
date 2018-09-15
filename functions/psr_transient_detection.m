function [onset,offset] = psr_transient_detection(signal,threshold)

% PSR_TRANSIENT_DETECTION - Detects signal transients
%
% Syntax:  [onset,offset] = psr_transient_detection(signal,threshold)
%
% Inputs:
%    signal    - Time series data
%    threshold - Threshold for detecting transient
%
% Outputs:
%    onset  -  Onsets for all transients
%    offset - Offsets for all transients

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Normalize
signal = signal / max(signal);
signal = single(signal); % convert to single type

% Find transients
signalDiff = diff(signal);

% On- and offsets
transitions = find(abs(signalDiff) > threshold);

% Join neighbouring transitions
d = diff(transitions);
transitions = transitions(find(d ~= 1) + 1);

% Check if we are dealing with onset or offset
amplitudes = signalDiff(transitions);
amplitudes = amplitudes ./ abs(amplitudes);

id = find(diff(amplitudes) == -2);
onset  = transitions(id);
offset = transitions(id + 1);

end