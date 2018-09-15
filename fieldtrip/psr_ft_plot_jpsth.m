function cfg = psr_ft_plot_jpsth(jpsth,parameters,unitIDs)

% PSR_FT_PLOT_JPSTH - Wrapper function for PSR_FT_SPIKE_PLOT_JPSTH
% This function calls the (modified) FieldTrip function:
% PSR_FT_SPIKE_PLOT_JPSTH, which plots the joint peri-stimulus time
% histogram (JPSTH)
%
% Syntax:  cfg = psr_ft_plot_jpsth(jpsth,parameters,unitIDs)
%
% Inputs:
%    jpsth      - Output from PSR_FT_JPSTH
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitIDs    - Unit IDs of the units we want to plot the JPSTH for
%
% Outputs:
%    cfg - Output from PSR_FT_SPIKE_PLOT_JPSTH
%
% See also: PSR_FT_SPIKE_PLOT_JPSTH

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Check input
if (isempty_field(jpsth,'jpsth.jpsth')); return; end

% Set input
cfg = [];
if (nargin == 3); cfg.channelcmb = psr_ft_combis(jpsth,unitIDs); end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.barplot'));     cfg.barplot     = parameters.analysis.jpsth.plot.barplot;     end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.colorbar'));    cfg.colorbar    = parameters.analysis.jpsth.plot.colorbar;    end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.colormap'));    cfg.colormap    = parameters.analysis.jpsth.plot.colormap;    end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.interpolate')); cfg.interpolate = parameters.analysis.jpsth.plot.interpolate; end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.window'));      cfg.window      = parameters.analysis.jpsth.plot.window;      end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.gaussvar'));    cfg.gaussvar    = parameters.analysis.jpsth.plot.gaussvar;    end
if (~isempty_field(parameters,'parameters.analysis.jpsth.plot.winlen'));      cfg.winlen      = parameters.analysis.jpsth.plot.winlen;      end

% Call the plotting function
cfg = psr_ft_spike_plot_jpsth(cfg,jpsth);

end