%% Stimuli

parameters.analysis.stimoffset = 0;

% Windows to extract before and after stimulus [s] 
parameters.analysis.spk.t_win = [-0.25,0.50]; % Spiking window
parameters.analysis.lfp.t_win = [-0.25,0.50]; % Local field potential window (should include some padding)

%% Fast fourier transform

%% psr_lfp_fft

parameters.analysis.fft.taper      = 'hanning';  % Taper window type ('dpss','hanning',...)
parameters.analysis.fft.pad        = 'nextpow2'; % 
parameters.analysis.fft.foi        = 0.5:0.1:70; % Frequencies of interest [Hz]
parameters.analysis.fft.tapsmofrq  = 2;          % Amount of spectral smoothing through multi-tapering (when taper = 'dpss') [Hz]
parameters.analysis.fft.keepchans  = false;      % Whether to keep individual channels, or to average over them

%% psr_lfp_plot_fft

parameters.analysis.fft.plot.alpha   = 0.2;
parameters.analysis.fft.plot.color   = 'k';       
parameters.analysis.fft.plot.error   = 'std';     % Plot error range: 'no','std' (default), 'sem'
parameters.analysis.fft.plot.powtype = 'decibel'; % 
parameters.analysis.fft.plot.smfreq  = [];        % Gaussian smoothing kernel standard deviation (no smoothing when <= 0) [Hz]

%% Time-Frequency analysis

%% psr_lfp_tfa

parameters.analysis.tfa.taper      = 'hanning';
parameters.analysis.tfa.pad        = 'nextpow2';
parameters.analysis.tfa.foi        = 0.5:0.5:70; % Frequencies of interest [Hz]
parameters.analysis.tfa.tapsmofrq  = 2;          % Amount of spectral smoothing through multi-tapering (when taper = 'dpss') [Hz]
parameters.analysis.tfa.toi        = 'all';      % The times on which the analysis windows should be centered [s]
parameters.analysis.tfa.ncycles    = 5;          % Number of cycles for frequency in sliding time window. If empty, "t_ftimwin" is used instead.
parameters.analysis.tfa.t_ftimwin  = 2;          % Length of sliding time window [s]
parameters.analysis.tfa.keepchans  = false;      % Whether to keep individual channels, or to average over them

%% psr_lfp_baseline

parameters.analysis.tfa.baseline     = [];        % 'pre' or [begin end] in seconds
parameters.analysis.tfa.baselineType = 'decibel'; % absolute', 'relative', 'relchange', 'normchange' or 'decibel'

%% psr_lfp_plot_tfa

parameters.analysis.tfa.plot.powtype  = 'decibel'; % 'absolute' or 'decibel'
parameters.analysis.tfa.plot.colormap = jet(256);

%% Bandpower calculation

%% psr_lfp_bpw

parameters.analysis.bpw.trange  = [0 0.5]; % Time window for bandpower calculation [s]
parameters.analysis.bpw.maxmiss = 0.05;    % Maximum missing data  
parameters.analysis.bpw.fnorm   = [min(parameters.analysis.tfa.foi) max(parameters.analysis.tfa.foi)];
parameters.analysis.bpw.frange  = [...
    1,   4; ... % Delta
    4,   8; ... % Theta
    8,  12; ... % Alpha
    13, 30; ... % Beta
    30, 70];    % Low gamma

%% psr_lfp_plot_bpw

parameters.analysis.bpw.plot.alpha = 0.2;
parameters.analysis.bpw.plot.color = 'k';
parameters.analysis.bpw.plot.error = 'std';

%% Spike analysis

%% psr_ft_isi
parameters.analysis.isi.bins  = 0:0.0005:0.1; % ISI time bins [ms]
parameters.analysis.isi.param = 'coeffvar'; % compute the coefficient of variation (sd/mn of isis)

%% psr_ft_plot_isih
parameters.analysis.isi.plot.scatter  = 'no';
parameters.analysis.isi.plot.colormap = hot(256);
parameters.analysis.isi.plot.window   = 'gausswin';

%% psr_ft_psth
parameters.analysis.psth.binsize    = 0.005;  % [s]
parameters.analysis.psth.outputunit = 'rate';
parameters.analysis.psth.crop       = true;   % Crop trials to same length

%% psr_ft_plot_psth

%% psr_ft_plot_raster

%% psr_ft_jpsth
parameters.analysis.jpsth.shuffle    = true;  % Shuffle trial correction
parameters.analysis.jpsth.normalize  = 'yes'; % Normalize JPSTH using method from van Aertsen et al. (1989)
parameters.analysis.jpsth.keeptrials = 'no';

%% psr_ft_plot_jpsth
parameters.analysis.jpsth.plot.colormap = jet(256);

%% psr_ft_spike_xcorr
parameters.analysis.xcorr.shuffle    = true;  % Shuffle trial correction
parameters.analysis.xcorr.maxlag     = 0.1;   % Maximum lag for correlogram [s]
parameters.analysis.xcorr.outputunit = 'proportion';
parameters.analysis.xcorr.binsize    = 0.001; % [s]
parameters.analysis.xcorr.debias     = 'yes';