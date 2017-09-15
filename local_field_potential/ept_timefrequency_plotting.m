function ept_timefrequency_plotting(freq,parameters,savePath,filename,baseline)

if (nargin < 3); savePath = []; end % save in current working directory
if (nargin < 4); filename = []; end
if (nargin < 5); baseline = []; end

filename = ['Freq' filename(7:end)];

cfg = [];
cfg.trials    = 'all';
cfg.maskstyle = 'saturation';
baselinetype  = 'relative';

if (isempty(baseline))
    cfg.baseline     = [parameters.lfp.base_onset parameters.lfp.base_offset];
    cfg.baselinetype = baselinetype;
else
    cfg.baseline = 'no'; % use baseline from control session
    data = freq.powspctrm;
    if     (strcmp(baselinetype, 'absolute'));   data =  data -  baseline;
    elseif (strcmp(baselinetype, 'relative'));   data =  data ./ baseline;
    elseif (strcmp(baselinetype, 'relchange'));  data = (data -  baseline) ./ baseline;
    elseif (strcmp(baselinetype, 'normchange')); data = (data -  baseline) ./ (data + baseline);
    elseif (strcmp(baselinetype, 'db'));         data = 10 * log10(data ./ baseline);
    else error('unsupported method for baseline normalization: %s', baselinetype);
    end
    freq.powspctrm = data;
end

figure; set(gcf,'position',get(0,'screensize'));
ft_singleplotTFR(cfg,freq);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
h = colorbar;
ylabel(h, 'Relative power')
export_fig([savePath filename]);

end