function h = psr_lfp_plot_bpw(bandpwr,parameters)

% PSR_LFP_PLOT_BPW - Plot bandpower figure
%
% Syntax:  h = psr_lfp_plot_bpw(bandpwr,parameters)
%
% Inputs:
%    bandpwr    - Bandpower array, output from PSR_LFP_BANDPOWER or PSR_LFP_BPW
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    h - Handle for bandpower plot
%
% See also: PSR_LFP_BANDPOWER, PSR_LFP_BPW

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

nTrials = size(bandpwr,1);

alpha  = parameters.analysis.bpw.plot.alpha;
color  = parameters.analysis.bpw.plot.color;
frange = parameters.analysis.bpw.frange;
nrange = size(frange,1);

x = 1:nrange;
y = mean(bandpwr,1);
z = [];

switch parameters.analysis.bpw.plot.error
    case 'std'; z = std(bandpwr,[],1);
    case 'sem'; z = std(bandpwr,[],1) / sqrt(nTrials);
end

x = x';
y = y';
z = z';

h = plot(x,y,'-o','Color',color,'MarkerEdgeColor',color,'MarkerFaceColor',color);
if ~isempty(z)
    hold on;
    plot_ci(x,[y,(y - z),(y + z)], ...
        'PatchColor',    color,    ...
        'PatchAlpha',    alpha,    ...
        'MainLineStyle', '-',      ...
        'MainLineWidth', 1.0,      ...
        'MainLineColor', color,    ...
        'LineColor',     color,    ...
        'LineStyle',     'none');
    hold off;
end

xlabel('$\bf{Frequency [Hz]}$',  'Interpreter','Latex');
ylabel('$\bf{Relative \ power}$','Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

xlim([min(x)-1 max(x)+1]);
xticks(x);

xTicks = cell(0,0);
for iRange = 1:nrange
    xTicks{1,iRange} = [num2str(frange(iRange,1)) ' - ' num2str(frange(iRange,2))];
end
xticklabels(xTicks);

if (~isempty_field(parameters,'parameters.analysis.bpw.plot.plim')); ylim(parameters.analysis.bpw.plot.plim); end

end