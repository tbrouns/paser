function psr_ft_plot_psth(psth,parameters,unitID)

% PSR_FT_PLOT_PSTH - Wrapper function for PSR_FT_SPIKE_PLOT_PSTH This
% This function calls the (modified) FieldTrip function:
% PSR_FT_SPIKE_PLOT_PSTH, which plots the peri-stimulus time histogram
% (PSTH)
%
% Syntax:  psr_ft_plot_psth(psth,parameters,unitID)
%
% Inputs:
%    psth       - Output from PSR_FT_PSTH
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitID     - Unit ID to plot the PSTH for

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

cfg = [];
if (nargin == 3); cfg.spikechannel = psth.label{unitID}; end
if (~isempty_field(parameters,'parameters.analysis.psth.plot.errorbars')); cfg.errorbars = parameters.analysis.psth.plot.errorbars; end
if (~isempty_field(parameters,'parameters.analysis.psth.plot.ylim'));      cfg.ylim      = parameters.analysis.psth.plot.ylim;      end

psr_ft_spike_plot_psth(cfg,psth);
                
end