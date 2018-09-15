% PSR_PARAMETERS_ANALYSIS - Script to set default analysis parameters
% Relevant function name is given in each box. See comments before and
% after each parameter for its purpose.
% 
% For certain function that relate to the FieldTrip toolbox, the FieldTrip
% documentation of the related FieldTrip function (given below the function
% name) can also be consulted for more info on the various input
% parameters. Each field has the same function as the equivalent field in
% the FieldTrip documentation. For example, the parameter
% "parameters.analysis.psth.binsize" used in PSR_FT_PSTH is the same
% parameter as "cfg.binsize" in FT_SPIKE_PSTH
% 
% Syntax:  psr_parameters_analysis
%
% Outputs:
%    parameters - Structure that contains parameters for data analysis
% 
% See also: PSR_PARAMETERS_LOAD

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

%% PSR_STIMULUS_WINDOW
parameters.analysis.stimoffset = 0; % Offset to stimulus timestamps [sec]
parameters.analysis.spk.t_win  = [-0.25,0.50]; % Spiking window to extract around stimulus

%% PSR_FT_CONVERT2TRIALS
parameters.analysis.lfp.t_win  = [-2.00,4.00]; % LFP window to extract around stimulus (should include sufficient padding to avoid FieldTrip warnings)

%% Fast fourier transform

%% PSR_LFP_FFT
% Related FieldTrip function: FT_FREQANALYSIS 
parameters.analysis.fft.taper      = 'hanning';  % Taper window type ('dpss','hanning',...)
parameters.analysis.fft.pad        = 'nextpow2'; % Method for determining the padding size
parameters.analysis.fft.foi        = 0.5:0.1:70; % Frequencies of interest [Hz]
parameters.analysis.fft.tapsmofrq  = 2;          % Amount of spectral smoothing through multi-tapering (when taper = 'dpss') [Hz]
parameters.analysis.fft.keepchans  = false;      % Whether to keep individual channels, or to average over them

%% PSR_LFP_PLOT_FFT
parameters.analysis.fft.plot.alpha   = 0.2;       % Plot alpha (transparancy)
parameters.analysis.fft.plot.color   = 'k';       % Plot color
parameters.analysis.fft.plot.error   = 'std';     % Plot error range: 'no','std' (default), 'sem'
parameters.analysis.fft.plot.powtype = 'decibel'; % Power scale (empty or 'decibel') 
parameters.analysis.fft.plot.smfreq  = [];        % Gaussian smoothing kernel standard deviation (no smoothing when <= 0) [Hz]

%% Time-Frequency analysis

%% PSR_LFP_TFA
% Related FieldTrip function: FT_FREQANALYSIS 
parameters.analysis.tfa.taper      = 'hanning';  % Taper window type ('dpss','hanning',...)
parameters.analysis.tfa.pad        = 'nextpow2'; % Method for determining the padding size
parameters.analysis.tfa.foi        = 0.5:0.5:70; % Frequencies of interest [Hz]
parameters.analysis.tfa.tapsmofrq  = 2;          % Amount of spectral smoothing through multi-tapering (when taper = 'dpss') [Hz]
parameters.analysis.tfa.toi        = 'all';      % The times on which the analysis windows should be centered [s] (Note that you typically DO NOT want to use the default option here, due to very long processing times)
parameters.analysis.tfa.ncycles    = 5;          % Number of cycles for frequency in sliding time window. If empty, "t_ftimwin" is used instead.
parameters.analysis.tfa.t_ftimwin  = 2;          % Length of sliding time window [s]
parameters.analysis.tfa.keepchans  = false;      % Whether to keep individual channels, or to average over them

%% PSR_LFP_BASELINE
parameters.analysis.tfa.base.window  = 'pre';     % The baseline window. Given as: [begin end] in seconds, or set to 'pre' for everything before t = 0.
parameters.analysis.tfa.base.type    = 'decibel'; % The baseline calculation method: 'absolute', 'relative', 'relchange', 'normchange' or 'decibel'
parameters.analysis.tfa.base.ncycles = [];        % Baseline window defined by the number of cycles, i.e. frequency dependent. Set to empty if you don't want to use it, otherwise needs to be an integer.
                                                  % For example, if we set ncycles = 5, we would use a baseline window of 0.5 sec for the 10 Hz band, but a baseline window of only 0.25 sec for 20 Hz. 
                                                  
%% PSR_LFP_PLOT_TFA
% Related FieldTrip function: FT_SINGLEPLOTTFR
parameters.analysis.tfa.plot.powtype  = 'decibel'; % Power scale to use: 'absolute' or 'decibel'
parameters.analysis.tfa.plot.colormap = jet(256);

%% Bandpower calculation

%% PSR_LFP_BPW
parameters.analysis.bpw.trange  = [0, 0.5]; % Time window for bandpower calculation [s]
parameters.analysis.bpw.maxmiss = 0.05;    % Maximum missing data  
parameters.analysis.bpw.fnorm   = [min(parameters.analysis.tfa.foi) max(parameters.analysis.tfa.foi)];
parameters.analysis.bpw.frange  = [...
    1,   4; ... % Delta
    4,   8; ... % Theta
    8,  12; ... % Alpha
    13, 30; ... % Beta
    30, 70];    % Low gamma

%% PSR_LFP_PLOT_BPW
parameters.analysis.bpw.plot.alpha = 0.2;
parameters.analysis.bpw.plot.color = 'k';
parameters.analysis.bpw.plot.error = 'std';

%% Spike analysis

%% PSR_FT_ISI
% Related FieldTrip function: FT_SPIKE_ISI 
parameters.analysis.isi.bins  = 0:0.0005:0.1; % ISI time bins [ms]
parameters.analysis.isi.param = 'coeffvar';   % Compute the coefficient of variation (sd/mn of isis)

%% PSR_FT_PLOT_ISIH
% Related (modified) FieldTrip function: PSR_FT_SPIKE_PLOT_ISIRETURN
parameters.analysis.isi.plot.scatter  = 'no';
parameters.analysis.isi.plot.colormap = hot(256);
parameters.analysis.isi.plot.window   = 'gausswin';

%% PSR_FT_PSTH
% Related FieldTrip function: FT_SPIKE_PSTH
parameters.analysis.psth.binsize    = 0.005;  % The bin size of each PSTH bin [s]
parameters.analysis.psth.outputunit = 'rate'; % Which PSTH type to calculate
parameters.analysis.psth.crop       = true;   % Crop trials to same length

%% PSR_FT_PLOT_PSTH
% Related (modified) FieldTrip function: PSR_FT_SPIKE_PLOT_PSTH

%% PSR_FT_PLOT_RASTER
% Related (modified) FieldTrip function: PSR_FT_SPIKE_PLOT_RASTER

%% PSR_FT_JPSTH
% Related (modified) FieldTrip function: PSR_FT_SPIKE_JPSTH
parameters.analysis.jpsth.shuffle    = true;  % Shuffle trial correction
parameters.analysis.jpsth.normalize  = 'yes'; % Normalize JPSTH using method from van Aertsen et al. (1989)
parameters.analysis.jpsth.keeptrials = 'no';

%% PSR_FT_PLOT_JPSTH
% Related (modified) FieldTrip function: PSR_FT_SPIKE_PLOT_JPSTH
parameters.analysis.jpsth.plot.colormap = jet(256);

%% PSR_FT_XCORR
% Related (modified) FieldTrip function: PSR_FT_SPIKE_XCORR
parameters.analysis.xcorr.shuffle    = true;  % Shuffle trial correction
parameters.analysis.xcorr.maxlag     = 0.1;   % Maximum lag for correlogram [s]
parameters.analysis.xcorr.outputunit = 'proportion';
parameters.analysis.xcorr.binsize    = 0.001; 
parameters.analysis.xcorr.debias     = 'yes';

%% PSR_MULTISCALE_RELEVANCE
parameters.analysis.msr.dt_min  = 0.001; % Minimum time scale 
parameters.analysis.msr.dt_step = 100;   % Number of time scale steps between the minimum and maximum time scale
