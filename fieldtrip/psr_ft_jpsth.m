function [jpsth,jpsthShuff] = psr_ft_jpsth(psth,parameters,unitIDs)

% PSR_FT_JPSTH - Wrapper function for FieldTrip's PSR_FT_SPIKE_JPSTH
% This function calls the (modified) FieldTrip function:
% PSR_FT_SPIKE_JPSTH, which calculates the joint peri-stimulus time
% histogram (JPSTH) data
%
% Syntax:  [jpsth,jpsthShuff] = psr_ft_jpsth(psth,parameters,unitIDs)
%
% Inputs:
%    psth       - Output from PSR_FT_PSTH
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitIDs    - (Optional) Unit IDs of the units we want to calculate the
%                 JPSTH for
%
% Outputs:
%    jpsth      - Output from PSR_FT_SPIKE_JPSTH
%    jpsthShuff - Output from PSR_FT_SPIKE_JPSTH (shift predictor)
%
% See also: PSR_FT_SPIKE_JPSTH

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

cfg        = [];
cfg.method = 'jpsth';

if (nargin == 3); cfg.channelcmb = psr_ft_combis(psth,unitIDs); end
if (~isempty_field(parameters,'parameters.analysis.jpsth.keeptrials')); cfg.keeptrials    = parameters.analysis.jpsth.keeptrials; end
if (~isempty_field(parameters,'parameters.analysis.jpsth.normalize'));  cfg.normalization = parameters.analysis.jpsth.normalize;  end

jpsth = psr_ft_spike_jpsth(cfg,psth);

% Subtract the shift predictor
jpsthShuff = [];
if (parameters.analysis.jpsth.shuffle)
    cfg.method = 'shiftpredictor';
    jpsthShuff  = psr_ft_spike_jpsth(cfg,psth);
    jpsth.jpsth = jpsth.jpsth - jpsthShuff.shiftpredictor;
end

end