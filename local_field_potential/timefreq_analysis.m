close all
clear all

% Filtering sometimes returns following error: "Calculated filter
% coefficients have poles on or outside the unit circle and will not be
% stable. Try a higher cutoff frequency or a different type/order of
% filter." Therefore we must try different filter orders and cut-off
% frequencies until one works. Sub-sampling data also works.

% tbase_min   = 1;   % start time window for baseline (sec)
% tbase_max   = 99;  %   end time window for baseline (sec)
% tmin        = 100; % start time window for analysis (sec)
% tmax        = 200; %   end time window for analysis (sec)
% tstep       = 0.3; %      time bin size (sec)
% fstep       = 1.0; % frequency bin size (Hz)

fig1 = figure;  set(gcf,'position',get(0,'screensize'));
fig2 = figure;  set(gcf,'position',get(0,'screensize'));

% Spectrum parameters

twin        = 1.0; % length of time window (sec)
fmin        = 1.0; % frequency range min (Hz)
fmax        = 60;  % frequency range max (Hz)

Fbands = [1 4 7 15 31 fmax];
nfreq  = length(Fbands) - 1;

% Filter parameters

N      = 3:12;          % Try different filter orders
F_low  = 0.1:0.1:fmin;  % Try different frequencies (in Hz)
F_high = 300;

resamplingFactor = 3;

ADD_NOISE = 0;

% Filter parameters

cfg            = [];
cfg.channel    = 'all';
cfg.continuous = 'yes';

cfg.preproc.bpfilter   = 'yes';

files       = dir('data_AR_RS_*');
files       = char(files.name);
numfiles    = length(files(:,1));

for iTetrode = 1:numfiles
    
    data_input = importdata(files(iTetrode,:));
    
    % Re-sampling
    
    data        = data_input.trial{1};
    timestamps  = data_input.time{1};
    
    F_resample = resamplingFactor * F_high;
    
    if (F_resample < data_input.fsample) % only resample if new frequency is lower
        
        [data,timestamps] = resample(data',timestamps, F_resample); % resample
        
        data_input.trial        = {data'};
        data_input.time         = {timestamps};
        data_input.fsample      = F_resample;
        data_input.sampleinfo   = [1 length(timestamps)];
        
    end
    
    % Plot raw re-sampled signal
    % figure;
    % plot(timestamps,data');
    
    %% Pre-processing
    
    SUCCESS = 0;
    itr_tot = length(N) * length(F_low);
    
    for i = 1:length(N);
        
        if (SUCCESS); break; end
        
        for j = 1:length(F_low)
            
            cfg.preproc.bpfreq     = [F_low(j), F_high];
            cfg.preproc.bpfiltord  = N(i);
            
            try
                data  = ft_preprocessing(cfg, data_input);
                SUCCESS = 1;
                break;
            catch ME
                itr = j + (i - 1) * length(F_low);
                disp([ME.message]);
                SUCCESS = 0;
            end
        end
    end
    
    %% Time-freq analysis
    
    if (SUCCESS)
        
        Fs              = data.fsample;
        interval        = Fs;
        overlap         = round(Fs * 0.95);
        nfft            = Fs;
        nsection        = 3 * Fs;
        nsection_tot    = data.sampleinfo(2);
        nsections       = floor(nsection_tot / nsection);
        
        iStart          = 1;
        
        for iChan = 1:4
            for iSection = 1:nsections
                dat       = data.trial{1,1}(iChan,iStart:iStart+nsection);
                [S,F,T,P] = spectrogram(dat-mean(dat),interval,overlap,nfft,Fs);
                                
                figure(fig1);
                imagesc(T,F,10*log10(P));
                colorbar
                axis xy
                ylim([0 70]);
                xlabel('Time [s]');
                ylabel('Frequency [Hz]');
                export_fig(['Spectrogram_' num2str(iTetrode) '_' num2str(iChan) '_' num2str(iSection)]);
                
                figure(fig2);
                clf
                ColorSet = varycolor(nfreq);
                hold on;
                for i = 1:nfreq
                    semilogy(T, mean(P(Fbands(i):Fbands(i+1), :)), 'Color',ColorSet(i,:),'LineWidth',2)
                end
                set(gca,'TickLabelInterpreter', 'latex');
                set(gca, 'YScale', 'log')
                xlabel('$$\bf{Time \ (s)}$$' ,          'Interpreter', 'Latex')
                ylabel('$$\bf{Power \ (\mu V^2)}$$' ,   'Interpreter', 'Latex')
                legend('Delta','Theta','Alpha','Beta','Gamma');
                export_fig(['Spectrum_' num2str(iTetrode) '_' num2str(iChan) '_' num2str(iSection) '_2']);

                iStart = iStart + nsection + 1;
            end
        end
        
        %         % http://www.fieldtriptoolbox.org/faq/how_can_i_do_time-frequency_analysis_on_continuous_data
        %
        %         % Plot filtered signal
        %         %     figure;
        %         %     plot(data.time{1},data.trial{1});
        %
        %         if (ADD_NOISE) % add noise
        %             fNoise = 25;   % Frequency [Hz]
        %             aNoise = 250;  % Amplitude
        %             noise  =  aNoise*sin(2*pi.*data.time{1}.*fNoise);
        %             noise  =  noise.*sin(2*pi.*data.time{1}.*0.005);
        %             noise  = repmat(noise,4,1);
        %             data.trial{1,1} = noise + data.trial{1,1};
        %         end
        %
        %         % Segment data in trials
        %         cfg.length      = twin;
        %         cfg.overlap     = 0;
        %         data_segmented  = ft_redefinetrial(cfg, data);
        %
        %         cfg.method     = 'mtmfft';
        %         cfg.taper      = 'hanning';
        %         cfg.pad        = 'nextpow2';
        %         cfg.foilim     = [fmin fmax];
        %         cfg.keeptrials = 'yes';
        %         freq_segmented = ft_freqanalysis(cfg, data_segmented);
        %
        %         begsample   = data_segmented.sampleinfo(:,1);
        %         endsample   = data_segmented.sampleinfo(:,2);
        %         time        = ((begsample+endsample)/2) / data_segmented.fsample;
        %
        %         freq_continuous           = freq_segmented;
        %         freq_continuous.powspctrm = permute(freq_segmented.powspctrm, [2, 3, 1]);
        %         freq_continuous.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
        %         freq_continuous.time      = time;             % add the description of the time dimension
        %
        %         cfg.baseline     = [time(1) time(end)];
        %         cfg.baselinetype = 'normchange';
        %
        %         figure;
        %         ft_singleplotTFR(cfg,freq_continuous);
        %
        %         set(gca,'TickLabelInterpreter',                     'latex');
        %         xlabel('$$\bf{Time \ (s)}$$'       ,'Interpreter',  'Latex')
        %         ylabel('$$\bf{Frequency \ (Hz)}$$' ,'Interpreter',  'Latex')
        %         c = colorbar;
        %         title(c,{'$\bf{Power \ change}$'} ,'Interpreter','Latex');
        %         set(c.Title, 'Units', 'Normalized', 'Position', [1.0, 1.02, 0]);
        %         title(['$Tetrode \ ' num2str(iTetrode) '$'],'Interpreter','Latex');
        %         pause(0.5);
        %
        %         %     % Time-frequency analysis parameters
        %         %
        %         %     cfg             = [];
        %         %     cfg.output      = 'pow';
        %         %     cfg.pad         = 'nextpow2';
        %         %     cfg.method      = 'mtmconvol';
        %         %     cfg.foi         = fmin:fstep:fmax;               % frequency range
        %         %     cfg.toi         = tmin:tstep:tmax;               % time range
        %         %
        %         %     switch cfg.method
        %         %         case 'mtmconvol'
        %         %             cfg.taper       = 'hanning';
        %         %             cfg.t_ftimwin   = ones(length(cfg.foi),1).*twin;
        %         %         case 'mtmfft'
        %         %             cfg.tapsmofrq   =  4; % Hz
        %         %             cfg.taper       = 'dpss';
        %         %         case 'wavelet'
        %         %             cfg.width       = 7;
        %         %     end
        %         %
        %         %     freq             = ft_freqanalysis(cfg, data);
        %         %
        %         %     cfg.baseline     = [tbase_min tbase_max];
        %         %     cfg.baselinetype = 'relchange';
        %         %     figure;
        %         %     ft_singleplotTFR(cfg, freq);
        %         %
        %         %     set(gca,'TickLabelInterpreter',                     'latex');
        %         %     xlabel('$$\bf{Time \ (s)}$$'       ,'Interpreter',  'Latex')
        %         %     ylabel('$$\bf{Frequency \ (Hz)}$$' ,'Interpreter',  'Latex')
        %         %     c = colorbar;
        %         %     title(c,{'$\bf{Relative \ power \ change}$'} ,'Interpreter','Latex');
        %         %     set(c.Title, 'Units', 'Normalized', 'Position', [1.0, 1.02, 0]);
        %         %
        %         %     nfreq    = length(freq.freq);
        %         %     ntime    = length(freq.time);
        %         %     pow_spec = zeros(nfreq,ntime);
        %         %
        %         %     for ii = 1:nfreq
        %         %         for jj = 1:ntime
        %         %             pow_spec(ii, jj) = freq.powspctrm(1, ii, jj);
        %         %         end
        %         %     end
        %         %
        %         %     ColorSet = varycolor(nfreq);
        %         %     figure
        %         %     hold on;
        %         %     for m=1:nfreq
        %         %         plot(cfg.toi, pow_spec(m, :), 'Color',ColorSet(m,:))
        %         %     end
        %         %
        %         %     set(gca,'TickLabelInterpreter', 'latex');
        %         %     xlabel('$$\bf{Time \ (s)}$$' ,          'Interpreter', 'Latex')
        %         %     ylabel('$$\bf{Power \ (\mu V^2)}$$' ,   'Interpreter', 'Latex')
        
    end
    
end
