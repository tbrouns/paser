%% Stimuli

parameters.analysis.stimOffset_1 = 0; % [ms]
parameters.analysis.stimOffset_2 = 0; % [ms]

parameters.analysis.t_win = [-250,500]; % Window to extract before and after stimulus [ms]
parameters.analysis.t_del = [-0.5,0.5]; % Window on both sides of stimulus to delete all spikes [ms]

%% Fast fourier transform

parameters.analysis.fft.taper      = 'hanning';
parameters.analysis.fft.pad        = 'nextpow2';
parameters.analysis.fft.citype     = 'std';
parameters.analysis.fft.smoothing  = false;
parameters.analysis.fft.freq.lower = 0.5; % Lower limit of frequency bins [Hz]
parameters.analysis.fft.freq.upper =  70; % Upper limit of frequency bins [Hz]
parameters.analysis.fft.freq.step  = 0.1; 
parameters.analysis.fft.keepchans  = false; % Whether to keep individual channels, or to average over them

%% Time-Frequency analysis

parameters.analysis.tfa.base_type  = 'absolute';
parameters.analysis.tfa.base_twin  = 'pre';
parameters.analysis.tfa.taper      = 'hanning';
parameters.analysis.tfa.pad        = 'nextpow2';
parameters.analysis.tfa.freq.lower = 0.5;    % Lower limit of frequency bins [Hz]
parameters.analysis.tfa.freq.upper = 70;     % Upper limit of frequency bins [Hz]
parameters.analysis.tfa.freq.step  = 0.5;    % Frequency bin size [Hz]
parameters.analysis.tfa.time.step  = 0.02;   % Time bin size [sec]
parameters.analysis.tfa.twin       = 2;      % Length of sliding time window [s]
parameters.analysis.tfa.ncycles    = 5;      % Number of cycles for frequency in sliding time window. If empty, 'twin' is used instead.
parameters.analysis.tfa.keepchans  = false;  % Whether to keep individual channels, or to average over them

% Multi-taper smoothing
parameters.analysis.tfa.smfreq  = 2;      % Amount of spectral smoothing through multi-tapering [Hz]
parameters.analysis.tfa.smratio = 0.5;    % Amount of spectral smoothing through multi-tapering as ratio of frequencies of interest
parameters.analysis.tfa.smtype  = 'freq'; % Smoothing with absolute frequency value ('freq') or ratio of frequencies of interest ('ratio')

%% Bandpower calculation
parameters.analysis.bpw.trange  = [0 0.5]; % Time window for bandpower calculation [s]
parameters.analysis.bpw.maxmiss = 0.05;    % Maximum missing data  
parameters.analysis.bpw.fnorm   = [parameters.analysis.tfa.freq.lower parameters.analysis.tfa.freq.upper];
parameters.analysis.bpw.frange = [...
    1,   4; ... % Delta
    4,   8; ... % Theta
    8,  12; ... % Alpha
    13, 30; ... % Beta
    30, 70];    % Low gamma

%% Spike analysis

parameters.analysis.cluster.minrate = 1; % Minimum firing rate

% Cluster classification
parameters.analysis.cluster.thresholds = true;
parameters.analysis.cluster.quality    = 3;

% Interspike interval histogram
parameters.analysis.isi.bins = 0:0.0005:0.1; % ISI time bins [ms]

% PSTH
parameters.analysis.psth.binsize = 0.005;

% JPSTH
parameters.analysis.jpsth.shuffle   = true; % Shuffle trial correction
parameters.analysis.jpsth.normalize = 'yes'; % Normalize JPSTH using method from van Aertsen et al. (1989)

% Correlogram
parameters.analysis.xcorr.shuffle = true; % Shuffle trial correction
parameters.analysis.xcorr.maxlag  = 0.1;  % Maximum lag for correlogram [s]