function XCorr = psr_ft_xcorr(spikesFT,parameters)

% PSR_FT_XCORR - Wrapper function for FieldTrip's PSR_FT_SPIKE_XCORR
% This function calls the FieldTrip function: PSR_FT_SPIKE_XCORR, which
% calculates cross-correlation metrics
%
% Syntax:  XCorr = psr_ft_xcorr(spikesFT,parameters)
%
% Inputs:
%    spikesFT   - FieldTrip spike structure
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    XCorr - Output from PSR_FT_SPIKE_XCORR
%
% See also: PSR_FT_SPIKE_XCORR, PSR_FT_CONVERT2FIELDTRIP

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

cfg            = [];
cfg.method     = 'xcorr';
cfg.keeptrials = 'yes';

if (~isempty_field(parameters,'parameters.analysis.xcorr.binsize'));    cfg.binsize    = parameters.analysis.xcorr.binsize;    end
if (~isempty_field(parameters,'parameters.analysis.xcorr.debias'));     cfg.debias     = parameters.analysis.xcorr.debias;     end
if (~isempty_field(parameters,'parameters.analysis.xcorr.maxlag'));     cfg.maxlag     = parameters.analysis.xcorr.maxlag;     end
if (~isempty_field(parameters,'parameters.analysis.xcorr.outputunit')); cfg.outputunit = parameters.analysis.xcorr.outputunit; end
    
XCorr = psr_ft_spike_xcorr(cfg,spikesFT);

% Subtract the shift predictor
if (parameters.analysis.xcorr.shuffle) 
    cfg.method  = 'shiftpredictor';
    Xshuff      = psr_ft_spike_xcorr(cfg,spikesFT);
    XCorr.xcorr = XCorr.xcorr - Xshuff.shiftpredictor;
end

end