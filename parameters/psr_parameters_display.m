%% DISPLAY PARAMETERS
parameters.display.metrics = true; % show metrics in plot title

% plot_waveforms
parameters.display.cmap           = hot(64); % Color map
parameters.display.waveform_ystep = 3;       % Steps in 2D waveform image given in units of signal  
parameters.display.waveform_clims = 0.35;    % Color-bar limit [Between: 0 - 1.0] 

% spike-train correlations parameters
parameters.display.show_isi                = 0;     % 0 = show autocorrelation by default, 1 = show ISI histogram
parameters.display.max_autocorr_to_display = 0.1;   % s, maximum time lag in autocorrelation/cross-correlation plots
parameters.display.max_isi_to_display      = 0.025; % s, maximum time lag to show in ISI histograms,
parameters.display.correlations_bin_size   = 2;     % ms, bin width used for spike time correlations
parameters.display.isi_bin_size            = 0.5;   % ms, bin width used for ISI histograms

% plot_stability parameters
parameters.display.stability_bin_size = 10;   % s, bin used for estimating firing rate
parameters.display.max_scatter        = 5000; % maximum number of waveforms to show in amplitude scatter plot
parameters.display.max_artifacts      = 50;   % longest artifact sections

% detection criterion
parameters.display.show_gaussfit = true;