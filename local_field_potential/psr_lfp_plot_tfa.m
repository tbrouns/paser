function h = psr_lfp_plot_tfa(data,parameters)

if (isempty(parameters.analysis.tfa.base.window))
    powtype = lower(parameters.analysis.tfa.plot.powtype);
    switch powtype; case 'decibel'; data.powspctrm = 10 * log10(data.powspctrm); end
else
    powtype = lower(parameters.analysis.tfa.base.type);
end
yLabelStr = [upper(powtype(1)) powtype(2:end) ' \ power'];

% Check data

nDims = ndims(data.powspctrm);
nTimes = length(data.time);
if (nTimes ~= size(data.powspctrm,nDims))
    if     (nDims == 3); data.powspctrm = data.powspctrm(:,:,  1:nTimes);
    elseif (nDims == 4); data.powspctrm = data.powspctrm(:,:,:,1:nTimes);
    end
end

% Spectrogram (Power vs. Frequency vs. Time)

% Remove artifacts field
data = psr_remove_field(data,'missing');
data = psr_remove_field(data,'artifacts');

if (all(isnan(data.powspctrm(:)))); return; end

% Plot spectrogram
cfg             = [];
cfg.colormap    = parameters.analysis.tfa.plot.colormap;
cfg.interactive = 'no';
cfg.fontsize    = 11; % Default font-size
ft_singleplotTFR(cfg,data); h = gca;
xlabel('$\bf{Time \ [s]}$',      'Interpreter','Latex');
ylabel('$\bf{Frequency \ [Hz]}$','Interpreter','Latex');

c = colorbar;
ylabel(c,['$\bf{' yLabelStr '}$'],'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

tmin = min(data.time);
tmax = max(data.time);
xlim([tmin tmax]);

title('');

cmap = colormap;
set(gca,'Color',cmap(round(0.5 * size(cmap,1)),:));

if (~isempty_field(parameters,'parameters.analysis.tfa.plot.tlim'));  xlim(parameters.analysis.tfa.plot.tlim); end
if (~isempty_field(parameters,'parameters.analysis.tfa.plot.flim'));  ylim(parameters.analysis.tfa.plot.flim); end
if (~isempty_field(parameters,'parameters.analysis.tfa.plot.plim')); caxis(parameters.analysis.tfa.plot.plim); end

end