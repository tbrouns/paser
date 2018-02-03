function h = psr_lfp_plotting(input)

% input.artifact:      Whether to plot artifact lines in spectrogram (true or false, default: false)
% input.artifactAlpha: Transparency level of artifact lines (between 0 and 1, default: 0.5)

% Grab data

data = input.freq;

% Initialize parameters

cfg = [];

if (psr_isempty_field(data,'time')); cfg.plotType = 'fft';
else,                                cfg.plotType = 'tfa';
end

cfg.artifact      = false;
cfg.artifactAlpha = 0.5;
cfg.baseline      = [];
cfg.ci            = false; % Plot error around power spectrum
cfg.citype        = 'std';
cfg.color         = 'k';
cfg.interactive   = 'no';
cfg.powtype       = 'absolute';
cfg.smoothing     = false; % Smoothing of power spectrum through Gaussian convolution

% Parse inputs

if (~psr_isempty_field(input,'artifact'));      cfg.artifact      = input.artifact;      end
if (~psr_isempty_field(input,'artifactAlpha')); cfg.artifactAlpha = input.artifactAlpha; end
if (~psr_isempty_field(input,'baseline'));      cfg.baseline      = input.baseline;      end
if (~psr_isempty_field(input,'powtype'));       cfg.powtype       = input.powtype;       end
if (~psr_isempty_field(input,'ci'));            cfg.ci            = input.ci;            end
if (~psr_isempty_field(input,'citype'));        cfg.citype        = input.citype;        end
if (~psr_isempty_field(input,'color'));         cfg.color         = input.color;         end
if (~psr_isempty_field(input,'plotType'));      cfg.plotType      = input.plotType;      end
if (~psr_isempty_field(input,'sigma'));         cfg.sigma         = input.sigma;         end
if (~psr_isempty_field(input,'smoothing'));     cfg.smoothing     = input.smoothing;     end
if (~psr_isempty_field(input,'xlim'));          cfg.xlim          = input.xlim;          end
if (~psr_isempty_field(input,'ylim'));          cfg.ylim          = input.ylim;          end
if (~psr_isempty_field(input,'zlim'));          cfg.zlim          = input.zlim;          end

% Baseline correction
cfg.powtype = lower(cfg.powtype);
if (~isempty(cfg.baseline))
    switch cfg.powtype
        case 'absolute';   data.powspctrm =  data.powspctrm -  cfg.baseline;
        case 'relative';   data.powspctrm =  data.powspctrm ./ cfg.baseline;
        case 'relchange';  data.powspctrm = (data.powspctrm -  cfg.baseline) ./ cfg.baseline;
        case 'normchange'; data.powspctrm = (data.powspctrm -  cfg.baseline) ./ (data.powspctrm + cfg.baseline);
        case 'decibel';    data.powspctrm = 10 * log10(data.powspctrm ./ cfg.baseline);
    end
else % Convert data if specified
    switch cfg.powtype
        case 'decibel'; data.powspctrm = 10 * log10(data.powspctrm);
    end
end

% Plot results

yLabelStr = [upper(cfg.powtype(1)) cfg.powtype(2:end) ' \ power'];

switch cfg.plotType
    
    case 'fft' % Power spectrum (Power vs. Frequency)
        
        powSpctrm = data.powspctrm;
        freqArray = data.freq;
        
        nDims   = ndims(powSpctrm);
        nTrials =  size(powSpctrm,1);
        
        if (nDims == 3) % Average over trials
            powSpctrmMean = mean(powSpctrm,1);
            powSpctrmMean = permute(powSpctrmMean,[2 3 1]);
            if (nTrials > 1 && cfg.ci)
                switch cfg.citype
                    case 'std'; n = 1;
                    case 'sem'; n = sqrt(nTrials);
                end
                sd = std(powSpctrm,[],1);
                sd = permute(sd,[2 3 1]);
                l = powSpctrmMean - sd / n;
                u = powSpctrmMean + sd / n;
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
        freqArray(del) = [];
        powSpctrm(del) = [];
        
        if (cfg.ci) % Plot interval
            l(del) = [];
            u(del) = [];
            plot_ci(freqArray,[l' u'], ...
                'PatchColor',    cfg.color,  ...
                'PatchAlpha',    cfg.artifactAlpha,  ...
                'LineStyle',     'none');
        end
        
        % Plot power specturm
        h = plot(freqArray,powSpctrm,'LineWidth',1.5,'Color',cfg.color);
        xlabel( '$\bf{Frequency \ [Hz]}$','Interpreter','Latex');
        ylabel(['$\bf{' yLabelStr '}$'],  'Interpreter','Latex');
        set(gca,'TickLabelInterpreter','Latex');
        
    case 'tfa' % Spectrogram (Power vs. Frequency vs. Time)
        
        % Remove artifacts field
        artifacts = data.artifacts;
        data = rmfield(data,'artifacts');
        
        cfg.baseline = []; % Don't use built-in baseline
        
        % Plot spectrogram
        ft_singleplotTFR(cfg,data);
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
        
        % Add artifact lines
        if (cfg.artifact)
            ylimits = ylim;
            nartifacts = size(artifacts,1);
            hold on;
            for i = 1:nartifacts
                F = area(artifacts(i,:),[1 1] .* ylimits(end),'EdgeColor','none','FaceColor','r');
                alpha(F, cfg.artifactAlpha)
            end
        end
end

if (~psr_isempty_field(cfg,'xlim')); xlim(cfg.xlim); end
if (~psr_isempty_field(cfg,'ylim')); ylim(cfg.ylim); end
if (~psr_isempty_field(cfg,'zlim')); zlim(cfg.zlim); end

end