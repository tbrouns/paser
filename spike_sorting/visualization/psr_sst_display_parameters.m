function spikes = psr_sst_display_parameters(spikes)

%% DISPLAY PARAMETERS

% plot_waveforms
spikes.params.display.default_waveformmode = 3;       % 1 = show all waveforms, 2 = show 95% waveform bands, 3 = show 2d histogram
spikes.params.display.time_scalebar        = 1;       % ms, length of scalebar in milliseconds on waveform plots
spikes.params.display.cmap                 = hot(64); % default color map
spikes.params.display.waveform_step        = 2;       % Steps in 2D waveform image given in units of signal    
spikes.params.display.default_outlier_method = 1; % 1 = estimate covariance from cluster, 2 = estimate covariance from background noise

% spike-train correlations parameters
spikes.params.display.show_isi                = 0;     % 0 = show autocorrelation by default, 1 = show ISI histogram
spikes.params.display.max_autocorr_to_display = 0.1;   % s, maximum time lag in autocorrelation/cross-correlation plots
spikes.params.display.max_isi_to_display      = 0.025; % s, maximum time lag to show in ISI histograms,
spikes.params.display.correlations_bin_size   = 2;     % ms, bin width used for spike time correlations
spikes.params.display.isi_bin_size            = 0.5;   % ms, bin width used for ISI histograms
spikes.params.display.default_xcorr_mode      = 1;     % 0 = show autocorrelation of merged spike trains by default, 1 = show x-correlation
spikes.params.display.trial_spacing           = 0.5;   % s, time padding between trials

% plot_stability parameters
spikes.params.display.stability_bin_size = 10;   % s, bin used for estimating firing rate
spikes.params.display.max_scatter        = 5000; % maximum number of waveforms to show in amplitude scatter plot

% figure layout
spikes.params.display.default_figure_size         = [0.05 0.1 0.9 0.8]; % location of default figure
spikes.params.display.figure_font_size            = 8;                  % default font size
spikes.params.display.initial_split_figure_panels = 4;

% axes/panels grid parameters
spikes.params.display.margin       = 140; % number of pixels between plots
spikes.params.display.outer_margin = 80;  % number of pixels around at figure margin
spikes.params.display.width        = 140; % number of pixels in plot width
spikes.params.display.aspect_ratio = 2/3; % aspect ratio of plots

end