function psr_sst_plot_poisson(spikes,clustID,parameters)

% Poisson distribution
id     = find([spikes.clusters.vars.id] == clustID);
twin   = parameters.cluster.stability_win;
frate  = spikes.clusters.vars(id).frate;
N      = spikes.clusters.vars(id).p_dist;
lambda = frate * twin;
x      = 0:length(N)-1;
y      = poisspdf(x,lambda);

% Plot

bar(x,N,1.0,'FaceColor','k','EdgeColor','none','FaceAlpha',1.0);
hold on;
plot(x,y,'r','LineWidth',1.5)
xlim([x(1) x(end)]);

ymax = 2 * max(y);
if (ymax > 1.0); ymax = 1.0; end
ylim([0 ymax]);

xlabelstr = ['\bf{No. \ of \ spikes \ in \ ' num2str(parameters.cluster.stability_win) '\ sec \ interval}'];
xlabel(   xlabelstr, 'Interpreter','Latex');
ylabel('\bf{Count}', 'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

ylabh = get(gca,'ylabel');
set(ylabh,'Units','normalized');
set(ylabh,'position', get(ylabh,'position') - [0.01 0.20 0]);

end
