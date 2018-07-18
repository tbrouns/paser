function output = psr_ft_plot_isih(isih,parameters,clustID,plotRaw)

cfg = [];

if (nargin < 3 || isempty(clustID)); clustID = 1; end % Can only plot one cluster at-a-time
if (nargin > 3 && plotRaw); parameters.analysis.isi.plot.window = 'no'; end

cfg.spikechannel = isih.label{clustID};
    
if (~isempty_field(parameters,'parameters.analysis.isi.plot.window'));   cfg.window   = parameters.analysis.isi.plot.window;   end 
if (~isempty_field(parameters,'parameters.analysis.isi.plot.winlen'));   cfg.winlen   = parameters.analysis.isi.plot.winlen;   end
if (~isempty_field(parameters,'parameters.analysis.isi.plot.colormap')); cfg.colormap = parameters.analysis.isi.plot.colormap; end 
if (~isempty_field(parameters,'parameters.analysis.isi.plot.scatter'));  cfg.scatter  = parameters.analysis.isi.plot.scatter;  end 

output = psr_ft_spike_plot_isireturn(cfg,isih);

end