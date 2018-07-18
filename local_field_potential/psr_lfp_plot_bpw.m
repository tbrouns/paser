function h = psr_lfp_plot_bpw(bandpwr,parameters)

nTrials = size(bandpwr,1);

alpha  = parameters.analysis.bpw.plot.alpha;
color  = parameters.analysis.bpw.plot.color;
frange = parameters.analysis.bpw.frange;

x = (1:size(frange,1));
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
xticklabels({...
    [num2str(frange(1,1)) ' - ' num2str(frange(1,2))],...
    [num2str(frange(2,1)) ' - ' num2str(frange(2,2))],...
    [num2str(frange(3,1)) ' - ' num2str(frange(3,2))],...
    [num2str(frange(4,1)) ' - ' num2str(frange(4,2))],...
    [num2str(frange(5,1)) ' - ' num2str(frange(5,2))]});

if (~isempty_field(parameters,'parameters.analysis.bpw.plot.plim')); ylim(parameters.analysis.bpw.plot.plim); end

end