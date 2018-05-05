function cfg = psr_ft_plot_jpsth(jpsth,parameters,clustIDs)

% Check input
if (psr_isempty_field(jpsth,'jpsth.jpsth')); return; end

% Set input
cfg = [];
if (nargin == 3); cfg.channelcmb = psr_ft_combis(jpsth,clustIDs); end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.barplot'));     cfg.barplot     = parameters.analysis.jpsth.plot.barplot;     end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.colorbar'));    cfg.colorbar    = parameters.analysis.jpsth.plot.colorbar;    end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.colormap'));    cfg.colormap    = parameters.analysis.jpsth.plot.colormap;    end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.interpolate')); cfg.interpolate = parameters.analysis.jpsth.plot.interpolate; end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.window'));      cfg.window      = parameters.analysis.jpsth.plot.window;      end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.gaussvar'));    cfg.gaussvar    = parameters.analysis.jpsth.plot.gaussvar;    end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.plot.winlen'));      cfg.winlen      = parameters.analysis.jpsth.plot.winlen;      end

% Call plot function
cfg = psr_ft_spike_plot_jpsth(cfg,jpsth);

end