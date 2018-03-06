function h = psr_lfp_plot_fft(input)

% Grab data

data = input.freq;

% Initialize parameters

cfg             = [];
cfg.alpha       = 0.2;
cfg.ci          = false; % Plot error around power spectrum
cfg.citype      = 'std';
cfg.color       = 'k';
cfg.interactive = 'no';
cfg.powtype     = 'absolute';
cfg.smoothing   = false; % Smoothing of power spectrum through Gaussian convolution

% Parse inputs

if (~psr_isempty_field(input,'input.alpha'));     cfg.alpha     = input.alpha;     end
if (~psr_isempty_field(input,'input.powtype'));   cfg.powtype   = input.powtype;   end
if (~psr_isempty_field(input,'input.ci'));        cfg.ci        = input.ci;        end
if (~psr_isempty_field(input,'input.citype'));    cfg.citype    = input.citype;    end
if (~psr_isempty_field(input,'input.color'));     cfg.color     = input.color;     end
if (~psr_isempty_field(input,'input.sigma'));     cfg.sigma     = input.sigma;     end
if (~psr_isempty_field(input,'input.smoothing')); cfg.smoothing = input.smoothing; end
if (~psr_isempty_field(input,'input.xlim'));      cfg.xlim      = input.xlim;      end
if (~psr_isempty_field(input,'input.ylim'));      cfg.ylim      = input.ylim;      end
if (~psr_isempty_field(input,'input.zlim'));      cfg.zlim      = input.zlim;      end

% Convert data if specified
cfg.powtype = lower(cfg.powtype);
switch cfg.powtype
    case 'decibel'; data.powspctrm = 10 * log10(data.powspctrm);
end

% Plot results

yLabelStr = [upper(cfg.powtype(1)) cfg.powtype(2:end) ' \ power'];

% Power spectrum (Power vs. Frequency)

powSpctrm = data.powspctrm;
freqArray = data.freq;

nDims   = ndims(powSpctrm);
nTrials =  size(powSpctrm,1);

if (nDims == 3) % Average over trials
    keep = ~isinf(sum(powSpctrm,3));
    powSpctrm = powSpctrm(keep(:),:,:);
    powSpctrmMean = nanmean(powSpctrm,1);
    powSpctrmMean = permute(powSpctrmMean,[2 3 1]);
    if (nTrials > 1 && cfg.ci)
        switch cfg.citype
            case 'std'; n = 1;
            case 'sem'; n = sqrt(nTrials);
        end
        sd = std(powSpctrm,[],1);
        sd = permute(sd,[2 3 1]);
        l = powSpctrmMean - (sd / n);
        u = powSpctrmMean + (sd / n);
    else
        cfg.ci = false;
    end
    powSpctrm = powSpctrmMean;
end

if (cfg.smoothing)
    df = mean(diff(freqArray));
    powSpctrm = psr_gauss_smoothing(powSpctrm,df,cfg.sigma);
    if (cfg.ci)
        l = psr_gauss_smoothing(l,df,cfg.sigma);
        u = psr_gauss_smoothing(u,df,cfg.sigma);
    end
end

% Remove NaNs and Infs

del = (isnan(powSpctrm) | isinf(powSpctrm));
if (all(del)); return; end
freqArray(del) = [];
powSpctrm(del) = [];

if (cfg.ci) % Plot interval
    l(del) = [];
    u(del) = [];
    plot_ci(freqArray,[l' u'],   ...
        'PatchColor', cfg.color, ...
        'PatchAlpha', cfg.alpha, ...
        'LineStyle',  'none');
end

% Plot power specturm
h = plot(freqArray,powSpctrm,'LineWidth',1.5,'Color',cfg.color);
xlabel( '$\bf{Frequency \ [Hz]}$','Interpreter','Latex');
ylabel(['$\bf{' yLabelStr '}$'],  'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex');

if (~psr_isempty_field(cfg,'cfg.xlim')); xlim(cfg.xlim); end
if (~psr_isempty_field(cfg,'cfg.ylim')); ylim(cfg.ylim); end
if (~psr_isempty_field(cfg,'cfg.zlim')); zlim(cfg.zlim); end

end