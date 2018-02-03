%% Stimuli

parameters.analysis.stimOffset_1 = 0; % [ms]
parameters.analysis.stimOffset_2 = 0; % [ms]

parameters.analysis.stimCombination = true;
parameters.analysis.stimArray = [...
    0,  0;   ...
    10, 40;  ...
    50, 80;  ...
    90, 120];

%% %%%% Local field potential analysis

parameters.analysis.lfp.citype = 'std';

%% PSR_LFP_TIMEFREQ

parameters.analysis.lfp.freqRange = [...
    1,   4; ... % Delta
    4,   8; ... % Theta
    8,  12; ... % Alpha
    13, 30; ... % Beta
    30, 70];    % Low gamma

parameters.analysis.lfp.window    = [-500,500];
parameters.analysis.lfp.base_win  = [-500,  0];
parameters.analysis.lfp.base_type = 'absolute';
parameters.analysis.lfp.base      = false;
parameters.analysis.lfp.plot_mean = false;

% (Time-)Frequency analysis

parameters.analysis.tfa.method = 'mtmconvol';
parameters.analysis.tfa.taper  = 'hanning';
parameters.analysis.tfa.pad    = 'nextpow2';

parameters.analysis.tfa.smfreq     = 2;      % Amount of spectral smoothing through multi-tapering [Hz]
parameters.analysis.tfa.smratio    = 0.5;    % Amount of spectral smoothing through multi-tapering as ratio of frequencies of interest
parameters.analysis.tfa.smtype     = 'freq'; % Smoothing with absolute frequency value ('freq') or ratio of frequencies of interest ('ratio')
parameters.analysis.tfa.freq.lower = 1.0;    % Lower limit of frequency bins [Hz]
parameters.analysis.tfa.freq.upper = 60;     % Upper limit of frequency bins [Hz]
parameters.analysis.tfa.freq.step  = 0.5;    % Frequency bin size [Hz]
parameters.analysis.tfa.time.step  = 0.02;   % Time bin size [sec]
parameters.analysis.tfa.twin       = 2;      % Length of sliding time window [s]
parameters.analysis.tfa.ncycles    = 5;      % Number of cycles for frequency in sliding time window. If empty, 'twin' is used instead.
parameters.analysis.tfa.keepchans  = true;   % Whether to keep individual channels, or to average over them

%% %%%% Spike analysis

%% Cluster classification

parameters.analysis.cluster.classifier = false;
parameters.analysis.cluster.thresholds = true;
parameters.analysis.cluster.quality    = 3;

%% Spikes

n = 10; m = 20;
parameters.analysis.t_bin   = 11; % [ms]
parameters.analysis.t_win   = [-parameters.analysis.t_bin * (0.5 + n), parameters.analysis.t_bin * (0.5 + m)]; % Window to extract before and after stimulus [ms]
parameters.analysis.t_del   = [-0.5,0.5]; % window on both sides of stimulus to delete all spikes [ms]
parameters.analysis.t_pad   = 50; % [ms]
parameters.analysis.t_array = ...
    [-500, -50;  ...
    0,      50;  ...
    50,    100;  ...
    100,   200;  ...
    200,   400]; % [ms]

parameters.analysis.Nspikes = 10000;

% Baseline

n = 8; m = 2;
parameters.analysis.base_win = [-parameters.analysis.t_bin * (0.5 + n),-parameters.analysis.t_bin * (0.5 + m)];

% PSTH
n = 5; m = 10;
parameters.analysis.psth_win = [-parameters.analysis.t_bin * (0.5 + n), parameters.analysis.t_bin * (0.5 + m)]; % [ms]
parameters.analysis.psth_bin = [0.5 * parameters.analysis.t_bin, 1.5 * parameters.analysis.t_bin];

% JPSTH

parameters.analysis.jpsth_win = [-100,400];

% Stimulus onset firing rate difference

parameters.analysis.diff_win = [50,100,150,200];

% ISI analysis

parameters.analysis.isi_bin = 1;   % ISI time bin [ms]
parameters.analysis.isi_max = 150; % Maximum ISI  [ms]

% ACF

parameters.analysis.acf_bin = 1;  % Bin size for spike binning [ms]
parameters.analysis.acf_max = 100; % [ms]
parameters.analysis.acf_win = 300; % [ms]

% XCorr

parameters.analysis.xcorr_max = 500; % [ms]
parameters.analysis.xcorr_bin = 1;   % [ms]

% Spike-triggered average

parameters.analysis.sta.win = 100; % [ms]

%

parameters.analysis.y_step_amp = 0.02; % max

parameters.analysis.smooth = false;
parameters.analysis.sigma  = 50; % smoothing [ms]

% Figures

parameters.analysis.fig_hght =  600;
parameters.analysis.fig_wdth = 1200;