function psr_lfp_plotting(freq,baseline,params)

if (nargin < 2 || isempty(baseline)); baseline = []; end

cfg = [];
cfg.baseline  = 'no';
cfg.trials    = 'all';
cfg.maskstyle = 'opacity';
baselinetype  = params.lfp.base_type;

if (~isempty(baseline))
    data = freq.powspctrm;
    if     (strcmp(baselinetype, 'absolute'));   data =  data -  baseline;
    elseif (strcmp(baselinetype, 'relative'));   data =  data ./ baseline;
    elseif (strcmp(baselinetype, 'relchange'));  data = (data -  baseline) ./ baseline;
    elseif (strcmp(baselinetype, 'normchange')); data = (data -  baseline) ./ (data + baseline);
    elseif (strcmp(baselinetype, 'db'));         data = 10 * log10(data ./ baseline);
    else, error('unsupported method for baseline normalization: %s', baselinetype);
    end
        
    freq.powspctrm = data;
end

nDims = ndims(freq.powspctrm);
if (nDims == 4)
    freq.powspctrm = nanmean(freq.powspctrm,1);
    freq.powspctrm = squeeze(freq.powspctrm);
end

ft_singleplotTFR(cfg,freq);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
h = colorbar;
ylabel(h, [upper(baselinetype(1)) baselinetype(2:end) ' power'])
xmin = min(freq.time);
xmax = max(freq.time);
xlim([xmin xmax]);

end