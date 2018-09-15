function psr_ft_plot_raster(spikesFT,psth,parameters,unitIDs)

% PSR_FT_PLOT_RASTER - Wrapper function for PSR_FT_SPIKE_PLOT_RASTER 
% This function calls the (modified) FieldTrip function:
% PSR_FT_SPIKE_PLOT_RASTER, which plots both a raster plot and
% peri-stimulus time histogram (PSTH)
%
% Syntax:  psr_ft_plot_raster(spikesFT,psth,parameters,unitIDs)
%
% Inputs:
%    spikesFT   - FieldTrip spike structure
%    psth       - Output from PSR_FT_PSTH
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitIDs    - Unit IDs to generate the raster plot for
%
% See also: PSR_FT_SPIKE_PLOT_RASTER

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

cfg = [];
cfg.interactive  = 'no';
if (~isempty(psth))
    if (nargin == 4); cfg.spikechannel = psth.label(unitIDs); end
    if (~isempty_field(parameters,'parameters.analysis.raster.topplotfunc')); cfg.topplotfunc = parameters.analysis.raster.topplotfunc; end
    if (~isempty_field(parameters,'parameters.analysis.raster.markersize'));  cfg.markersize  = parameters.analysis.raster.markersize;  end
    if (~isempty_field(parameters,'parameters.analysis.raster.errorbars'));   cfg.errorbars   = parameters.analysis.raster.errorbars;   end

    psr_ft_spike_plot_raster(cfg,spikesFT,psth);
end
end