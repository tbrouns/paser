function h = psr_lfp_plot_tfa(input)

% Grab data

data = input.freq;

% Initialize parameters

cfg             = [];
cfg.baseline    = false;
cfg.convert     = false;
cfg.interactive = 'no';
cfg.powtype     = 'absolute';

% Parse inputs

if (~psr_isempty_field(input,'input.baseline')); cfg.baseline  = input.baseline; end
if (~psr_isempty_field(input,'input.base_t0'));  cfg.base_t0   = input.base_t0;  end
if (~psr_isempty_field(input,'input.base_t1'));  cfg.base_t1   = input.base_t1;  end
if (~psr_isempty_field(input,'input.convert'));  cfg.convert   = input.convert;  end
if (~psr_isempty_field(input,'input.powtype'));  cfg.powtype   = input.powtype;  end
if (~psr_isempty_field(input,'input.xlim'));     cfg.xlim      = input.xlim;     end
if (~psr_isempty_field(input,'input.ylim'));     cfg.ylim      = input.ylim;     end
if (~psr_isempty_field(input,'input.zlim'));     cfg.zlim      = input.zlim;     end

% Baseline correction
cfg.powtype = lower(cfg.powtype);
if (cfg.baseline)
    data = psr_lfp_baseline(data,cfg);
end

if (cfg.convert)
    switch cfg.powtype
        case 'decibel'; data.powspctrm = 10 * log10(data.powspctrm);
    end
end

% Spectrogram (Power vs. Frequency vs. Time)
yLabelStr = [upper(cfg.powtype(1)) cfg.powtype(2:end) ' \ power'];

% Remove artifacts field
data = psr_remove_field(data,'missing');
data = psr_remove_field(data,'artifacts');

cfg.baseline = []; % Don't use built-in baseline

if (all(isnan(data.powspctrm(:)))); return; end 

% Plot spectrogram
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

if (~psr_isempty_field(cfg,'cfg.xlim')); xlim(cfg.xlim); end
if (~psr_isempty_field(cfg,'cfg.ylim')); ylim(cfg.ylim); end
if (~psr_isempty_field(cfg,'cfg.zlim')); zlim(cfg.zlim); end

end