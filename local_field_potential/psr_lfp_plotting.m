function psr_lfp_plotting(freq,cfg)

% Parse inputs

if (nargin < 2); cfg = []; end

if (isfield(cfg,'baseline')); baseline = cfg.baseline;
else,                         baseline = [];
end

if (isfield(cfg,'basetype')); basetype = cfg.basetype;
else,                         basetype = 'absolute';
end

if (isfield(cfg,'trials')); trials = cfg.trials;
else,                       trials = [];
end

if (isfield(cfg,'maskstyle')); maskstyle = cfg.maskstyle;
else,                          maskstyle = [];
end

if (isfield(cfg,'xlim')); xlimits = cfg.xlim;
else,                     xlimits = [];
end

if (isfield(cfg,'ylim')); ylimits = cfg.ylim;
else,                     ylimits = [];
end

if (isfield(cfg,'zlim')); zlimits = cfg.zlim;
else,                     zlimits = [];
end

if (~isempty(baseline))
    data = freq.powspctrm;
    if     (strcmpi(basetype, 'absolute'));   data =  data -  baseline;
    elseif (strcmpi(basetype, 'relative'));   data =  data ./ baseline;
    elseif (strcmpi(basetype, 'relchange'));  data = (data -  baseline) ./ baseline;
    elseif (strcmpi(basetype, 'normchange')); data = (data -  baseline) ./ (data + baseline);
    elseif (strcmpi(basetype, 'db'));         data = 10 * log10(data ./ baseline);
    else, error('unsupported method for baseline normalization: %s', basetype);
    end
        
    freq.powspctrm = data;
end

cfg = [];
cfg.baseline  = 'no';
if (~isempty(trials));    cfg.trials    = trials;    end
if (~isempty(maskstyle)); cfg.maskstyle = maskstyle; end
if (~isempty(xlimits));   cfg.xlim      = xlimits;   end
if (~isempty(ylimits));   cfg.ylim      = ylimits;   end
if (~isempty(zlimits));   cfg.zlim      = zlimits;   end

ft_singleplotTFR(cfg,freq);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
h = colorbar;
ylabel(h, [upper(basetype(1)) basetype(2:end) ' power'])
xmin = min(freq.time);
xmax = max(freq.time);
xlim([xmin xmax]);

end