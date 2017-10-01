function parameters = ept_parameter_config()

% General

parameters.general.ext         = '.continuous'; % Extension of data files
parameters.general.pattern     = 'CH'; % Pattern to look for in data files
parameters.general.nelectrodes = 4; % Number of electrodes per polytrode

% Hardware specifications

parameters.process.spikes = false;
parameters.process.lfp    = false;
parameters.process.mfa    = false;

% Spike detection

parameters.spikes.ref_period = 1.5; % ms, refractory period (for calculation refractory period violations)

% spike detection parameters

% method:
% 'auto'    : threshold calculated from background noise
% 'manual'  : user defined threshold 
% 'mad'     : median absolute deviation (default)

parameters.spikes.method = 'mad'; 

% thresh:
% If 'auto'   : number of standard deviations above background noise
% If 'manual' : absolute threshold in microvolts
% If 'mad'    : number of standard deviations above background noise

parameters.spikes.thresh      = 3.0;
parameters.spikes.window_size = 1.5; % ms, width of a spike
parameters.spikes.shadow      = 0.5; % ms, enforced dead region after each spike
parameters.spikes.cross_time  = 0.6; % ms, alignment point for peak of waveform
parameters.spikes.max_jitter  = 0.6; % ms, width of window used to detect peak after threshold crossing

% General artifact parameters
parameters.spikes.artifacts_corr   = 0.5;  % correlation threshold
parameters.spikes.artifacts_offset = 0.15; 

% Filtering
parameters.spikes.bp_high  = 6000;
parameters.spikes.bp_low   = 600;
parameters.spikes.bp_order = 10;

parameters.spikes.tsection = 60; % time of each individual section that is processed (in minutes)

parameters.spikes.artifacts_removal = 1; % remove artifacts? (default: true)
parameters.spikes.artifacts_combine = 0; % Use both ADC and CONTINUOUS files to detect artifacts (currently NOT recommended)

%% SPIKE SORTING

% Dictionary learning: Focused Mixture Model (FMM)

parameters.sorting.fmm.p     = 1e-4;  % Changes how aggresively to cluster, range 0-1 (0: less clustering, 1: more clustering)
parameters.sorting.fmm.k     = 5; 
parameters.sorting.fmm.kmax  = 16;    % Maximum number of clusters
parameters.sorting.fmm.align = false; % Whether to align the waveforms 

% Continuous Basis Pursuit (CBP)

parameters.sorting.cbp.path = 'E:\Google Drive\Lab notebook Terry Brouns\Software\MATLAB\spike_sorting\CBPSpikesortDemo-master';

% Mixture of drifting t-distribution model (MDT)

parameters.sorting.mdt.dims              = 12; % Number of PC dimensions                               (default: 12)
parameters.sorting.mdt.nu                = 7;  % t-distribution nu parameter (smaller = heavier tails) (default: 7)
parameters.sorting.mdt.q_perhour         = 2;  % Drift regularization (smaller = more smoothing)       (default: 2)  
parameters.sorting.mdt.timeframe_minutes = 1;  % Time frame duration (mostly a computational thing)    (default: 1)

% UltraMegaSort2000 (UMS)

parameters.sorting.ums.kmeans_size = 0.01; % target size for miniclusters as fraction of total number of spikes
parameters.sorting.ums.agg_cutoff  = 0.0001;  % higher = less aggregation, lower = more aggregation

% Super paramagnetic clustering 

parameters.sorting.spc.dims = 3;

%% CLUSTER PARAMETERS & QUALITY CONTROL

parameters.cluster.method  = 'fmm';

parameters.cluster.upper_rpv  = 0.10; % Maximum fraction of refractory period violations (RPVs) in cluster
parameters.cluster.upper_miss = 0.10; % Maximum fraction of missing spikes in cluster due to threshold

parameters.cluster.amplitude_max = 300; % maximum absolute mean spike amplitude of cluster (microvolt)
parameters.cluster.thresh        = 4;   % number of standard deviations above background noise for mean amplitude of highest amplitude channel
parameters.cluster.lag_max       = 5;   % how close maximum cross correlation needs to be to zero lag point
parameters.cluster.thresh_xcorr  = 1;   % threshold for channel amplitude in order to be considered for cross-correlation
parameters.cluster.size_min      = 10;  % minimum number of spikes in cluster (change to mean firing rate threshold)

parameters.cluster.outlier_chi   = 0.001; 
parameters.cluster.outlier_std   = 3;     % max num. of SDs from centre of PCA cluster

%% LOCAL FIELD POTENTIAL

parameters.lfp.rsfactor = 4; % Multiple of high band-pass filter frequency to resample at

parameters.lfp.bp_high  = 300;
parameters.lfp.bp_low   = 0.1;
parameters.lfp.bp_order =   5;

parameters.lfp.trial_onset  = -0.2; % sec
parameters.lfp.trial_offset =  0.5; % sec
parameters.lfp.trial_padding = 1.0; % time window before and after trial on- and offset (sec)

% Time frequency analysis

% see FT_FREQANALYSIS help

parameters.lfp.method = 'mtmconvol';
parameters.lfp.taper  = 'hanning';
parameters.lfp.pad    = 'nextpow2';

parameters.lfp.freq_lower = 2.0; % Hz
parameters.lfp.freq_upper = 60;  % Hz
parameters.lfp.freq_step  = 0.5; % Hz

parameters.lfp.time_step = 0.02; % sec

parameters.lfp.ncycles = 5; % number of cycles for frequency in sliding time window

parameters.lfp.base_onset  = -0.3; % sec
parameters.lfp.base_offset = -0.1; % sec

parameters.lfp.base_type = 'relative';

parameters.lfp.artifact_freq     = 30; % maximum artifact frequency (Hz)
parameters.lfp.artifact_tsection = 10; % window length in which we detect background noise (sec)
parameters.lfp.artifact_thresh   =  7; % number of standard deviations above background noise
parameters.lfp.artifact_frac     = .1; % fraction of STD for determining artifact on- and offset based on derivative 
parameters.lfp.artifact_padding  =  1; % in secs

% Magnetic field artifact (MFA) parameters

% parameters.mfa.thresh   = 4.0;  % stds above noise
% parameters.mfa.fmm_p    = 0.01; % MFA dictionary learning parameter
% parameters.mfa.freq_off = 0.02; % offset for magnetic field artifact periodicity (fraction of period)
% parameters.mfa.freq_max = 1.00; % Upper bound on MFA frequency (Hz) 
% parameters.mfa.freq_min = 0.50; % Lower bound on MFA frequency (Hz) 
% parameters.mfa.std_corr = 12.0;
% parameters.mfa.thr_corr = 0.05;
% parameters.mfa.clus_min = 0.25;
% parameters.mfa.range    =    5;
% parameters.mfa.frac     = 0.75;

parameters.mfa.p2ptime  = 0.25; % time between min and max peaks (ms)
parameters.mfa.bp_high  = 6000; % Hz
parameters.mfa.bp_low   = 1000; % Hz
parameters.mfa.bp_order = 10;
parameters.mfa.thresh   = 8; % # of STDs

parameters.mfa.control = 200; % number of windows taken from control data

end