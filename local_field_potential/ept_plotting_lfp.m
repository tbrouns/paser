function ept_plotting_lfp(freq,parameters,savePath,filename,baseline)

if (nargin < 3); savePath = []; end % save in current working directory
if (nargin < 4); filename = []; end
if (nargin < 5); baseline = []; end

filename = ['Freq' filename(7:end)];

cfg = [];
cfg.baseline  = 'no';
cfg.trials    = 'all';
cfg.maskstyle = 'opacity';
baselinetype  = parameters.lfp.base_type;

%     cfg.baseline     = [parameters.lfp.base_onset parameters.lfp.base_offset];
%     cfg.baselinetype = baselinetype;

if (~isempty(baseline))
    % From FieldTrip:
    data = freq.powspctrm;
    if     (strcmp(baselinetype, 'absolute'));   data =  data -  baseline;
    elseif (strcmp(baselinetype, 'relative'));   data =  data ./ baseline;
    elseif (strcmp(baselinetype, 'relchange'));  data = (data -  baseline) ./ baseline;
    elseif (strcmp(baselinetype, 'normchange')); data = (data -  baseline) ./ (data + baseline);
    elseif (strcmp(baselinetype, 'db'));         data = 10 * log10(data ./ baseline);
    else error('unsupported method for baseline normalization: %s', baselinetype);
    end
    freq.powspctrm = data;
    
    if (strcmp(baselinetype,'relative'))
        cfg.zlim = [0 5];
    end
else % TEMP
    pow = freq.powspctrm;
    pow = nanmean(pow,3);
    pow = repmat(pow,1,1,size(freq.time,2));
    freq.powspctrm = pow;
end

figure; set(gcf,'position',get(0,'screensize'));
ft_singleplotTFR(cfg,freq);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
h = colorbar;
ylabel(h, 'Relative power')
xmin = parameters.lfp.trial_onset;
xmax = parameters.lfp.trial_offset;
xlim([xmin xmax]);
export_fig([savePath filename]);


end