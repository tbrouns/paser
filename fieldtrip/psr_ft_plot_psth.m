function psr_ft_plot_psth(psth,parameters,clustID)

cfg = [];
if (nargin == 3); cfg.spikechannel = psth.label{clustID}; end
if (~isempty_field(parameters,'parameters.analysis.psth.plot.errorbars')); cfg.errorbars = parameters.analysis.psth.plot.errorbars; end
if (~isempty_field(parameters,'parameters.analysis.psth.plot.ylim'));      cfg.ylim      = parameters.analysis.psth.plot.ylim;      end

psr_ft_spike_plot_psth(cfg,psth);
                
end