function psr_ft_plot_xcorr(XCorr,parameters,unitIDs)

% PSR_FT_PLOT_XCORR - Plots correlogram
%
% Syntax:  psr_ft_plot_xcorr(XCorr,parameters,unitIDs)
%
% Inputs:
%    XCorr      - Output from PSR_FT_XCORR
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    unitIDs    - Two-element vector containing the unit IDs we want to
%                 calculate the correlogram for
%
% See also: PSR_FT_XCORR

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (length(unitIDs) ~= 2)
    str = {'ERROR in "psr_ft_plot_xcorr": Input must contain two unit indices.'};
    psr_show_warning(str);
    return;
end

% Initialize

cfg            = [];
cfg.color      = 'k';
cfg.linewidth  = 1.5;
cfg.marker     = 'none';
cfg.markersize = 3;

if (~isempty_field(parameters,'parameters.analysis.xcorr.plot.color'));      cfg.color      = parameters.analysis.xcorr.plot.color;      end
if (~isempty_field(parameters,'parameters.analysis.xcorr.plot.linewidth'));  cfg.linewidth  = parameters.analysis.xcorr.plot.linewidth;  end
if (~isempty_field(parameters,'parameters.analysis.xcorr.plot.marker'));     cfg.marker     = parameters.analysis.xcorr.plot.marker;     end
if (~isempty_field(parameters,'parameters.analysis.xcorr.plot.markersize')); cfg.markersize = parameters.analysis.xcorr.plot.markersize; end

t = XCorr.time;
y = squeeze(XCorr.xcorr(unitIDs(1),unitIDs(2),:));
plot(t,y,...
    'Color',     cfg.color,...
    'LineWidth', cfg.linewidth,...
    'Marker',    cfg.marker,...
    'MarkerSize',cfg.markersize);
   
ylabelStr = [];
switch XCorr.cfg.outputunit
    case 'raw';        ylabelStr = 'Cross-correlation';
    case 'center';     ylabelStr = 'Scaled \ cross-correlation';
    case 'proportion'; ylabelStr = 'Proportion \ of \ occurence';        
end
ylabel(['\bf{' ylabelStr '}'],'Interpreter','Latex');
xlabel( '\bf{Lag \ [s]}',     'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

end