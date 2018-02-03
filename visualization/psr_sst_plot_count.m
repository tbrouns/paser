function psr_sst_plot_count(spikes,clustID,parameters)

% Poisson distribution
I      = find([spikes.clusters.metrics.id] == clustID);
twin   = parameters.cluster.stability.win;
frate  = spikes.clusters.metrics(I).frate;
xDist  = spikes.clusters.metrics(I).cxDist;
yDist  = spikes.clusters.metrics(I).cyDist;
lambda = frate * twin;

if (isempty(xDist)); set(gca,'Visible','off'); return; end % Don't show plot
y = psr_sst_normpdf_stability(xDist,lambda,parameters);

% Plot

hold on;
bar(xDist,yDist,1.0,'FaceColor','k','EdgeColor','none','FaceAlpha',1.0);
plot(xDist,y,'r','LineWidth',1.5)

% Calculate plot limits

I = find(cumsum(yDist) > 0.99,1);
if (~isempty(I)); xmax = xDist(I);
else,             xmax = xDist(end);
end
ymax = 2 * max(y); if (ymax > 1.0); ymax = 1.0; end
xlim([-0.5, xmax + 0.5]);
ylim([0 ymax]);

% Set x-ticks

dx = mean(diff(get(gca,'xTick')));
if (dx < 1); xticks(0:xmax); end

xlabelstr = ['\bf{No. \ of \ spikes \ in \ ' num2str(parameters.cluster.stability.win) '\ sec \ interval}'];
xlabel(   xlabelstr, 'Interpreter','Latex');
ylabel('\bf{Count}', 'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

ylabh = get(gca,'ylabel');
set(ylabh,'Units','normalized');
set(ylabh,'position', get(ylabh,'position') - [0.01 0.20 0]);

end
