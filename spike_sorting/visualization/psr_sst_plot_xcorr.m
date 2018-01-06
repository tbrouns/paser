function psr_sst_plot_xcorr(spikes,clustID)

Fs = spikes.Fs;
id = logical([spikes.clusters.metrics.id] == clustID);
xc = spikes.clusters.metrics(id).xc;

if ~isempty(xc)
    n = 0.5 * (length(xc) - 1);
    lags = -n:n;
    lags = 1000 * lags / Fs;
    plot(lags,xc,'Color','k','LineWidth',2.0);
    xlabel('\bf{Time \ lag \ [ms]}',   'Interpreter','Latex');
    ylabel('\bf{Cross \ correlation}', 'Interpreter','Latex');
    set(gca,'TickLabelInterpreter','Latex');
    xlim = [min(lags) max(lags)];
    ylim([-1 1]);
    line([0 0], ylim,'Color','k','LineStyle',':');
    line( xlim,[0 0],'Color','k','LineStyle',':');
    axis([xlim ylim]);
else
    set(gca,'Visible','off'); 
end

end