function output = psr_ft_plot_isih(isih,parameters,unitID,RAW_PLOT)

% FUNCTION_NAME - Wrapper function for PSR_FT_SPIKE_PLOT_ISIRETURN
% This function calls the (modified) FieldTrip function:
% PSR_FT_SPIKE_PLOT_ISIRETURN, which plots the return plot (Poincare plot)
% for the inter-spike intervals (ISIs)
%
% Syntax:  output = psr_ft_plot_isih(isih,parameters,unitID,RAW_PLOT)
%
% Inputs:
%    isih       - Output from PSR_FT_ISI
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitID     - Unit ID to plot the return plot for
%    RAW_PLOT   - Boolean indicating whether to use smoothing, if:
%                 RAW_PLOT = true
%                 Then we use no smoothing
% Outputs:
%    output - Output from PSR_FT_SPIKE_PLOT_ISIRETURN
%
% See also: PSR_FT_SPIKE_PLOT_ISIRETURN

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

cfg = [];

if (nargin < 3 || isempty(unitID)); unitID = 1; end % Can only plot one cluster at-a-time
if (nargin > 3 && RAW_PLOT); parameters.analysis.isi.plot.window = 'no'; end

cfg.spikechannel = isih.label{unitID};
    
if (~isempty_field(parameters,'parameters.analysis.isi.plot.window'));   cfg.window   = parameters.analysis.isi.plot.window;   end 
if (~isempty_field(parameters,'parameters.analysis.isi.plot.winlen'));   cfg.winlen   = parameters.analysis.isi.plot.winlen;   end
if (~isempty_field(parameters,'parameters.analysis.isi.plot.colormap')); cfg.colormap = parameters.analysis.isi.plot.colormap; end 
if (~isempty_field(parameters,'parameters.analysis.isi.plot.scatter'));  cfg.scatter  = parameters.analysis.isi.plot.scatter;  end 

output = psr_ft_spike_plot_isireturn(cfg,isih);

end