function parameters = ept_parameter_config()

% General

parameters.general.ext         = '.continuous'; % Extension of data files
parameters.general.pattern     = 'CH'; % Pattern to look for in data files
parameters.general.nelectrodes = 4; % Number of electrodes per polytrode

% Pipeline parameters

parameters.process.spikes = 1; % Do spike detection
parameters.process.lfp    = 1; % Do LFP   detection
parameters.process.mfa    = 1; % Do MFA   detection

% Filtering

parameters.filter.fft_freq   = 10;   % Hz
parameters.filter.fft_pad    = 0.01; % Hz
parameters.filter.fft_thresh = 5;    % Standard deviations

% Spike detection

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
parameters.spikes.ref_period  = 1.5; % ms, refractory period (for calculation refractory period violations)

parameters.spikes.window_size = 1.5; % ms, width of a spike. Take samples symmetrically around peak
parameters.spikes.max_desync  = 0.2; % ms, maximum temporal desynchronization between channels of probe

% General artifact parameters
parameters.spikes.artifacts_corr   = 0.5;  % correlation threshold
parameters.spikes.artifacts_offset = parameters.spikes.max_desync; 

% Filtering
parameters.spikes.bp_high  = 6000;
parameters.spikes.bp_low   = 600;
parameters.spikes.bp_order = 10;

parameters.spikes.tsection = 60; % time of each individual section that is processed (in minutes)

parameters.spikes.artifacts_removal  = 1; % remove artifacts? (default: true)
parameters.spikes.artifacts_combine  = 0; % Use both ADC and CONTINUOUS files to detect artifacts (currently NOT recommended)
parameters.spikes.artifacts_subtract = 0; % subtract global median to remove artifacts

%% SPIKE SORTING

parameters.sorting.method  = 'KST';

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

% KiloSort (KST)

% Path directory should contain config file, named 'kilsort_config.m'
parameters.sorting.kst.path = 'E:\Google Drive\Lab notebook Terry Brouns\Software\MATLAB\spike_sorting\Repos\kilosort';
parameters.sorting.kst.cuda = 1; % use CUDA (GPU computing)

% Super paramagnetic clustering (SPC)

parameters.sorting.spc.dims = 3;

%% CLUSTER PARAMETERS & QUALITY CONTROL

% Thresholds

parameters.cluster.max_rpv       = 0.05; % Maximum fraction of refractory period violations (RPVs) in cluster
parameters.cluster.max_missing   = 0.05; % Maximum fraction of missing spikes in cluster due to threshold
parameters.cluster.max_amplitude = 300;  % Maximum absolute mean spike amplitude of cluster (microvolt)
parameters.cluster.max_artifact  = 1.25; % Maximum ratio between actual and expected number of spikes in LFP artifact region
parameters.cluster.max_lratio    = 1.0;  % Maximum L-ratio
parameters.cluster.max_lag       = parameters.spikes.max_desync; % Maximum lag of maximum pairwise cross correlation (ms)

parameters.cluster.min_isodist   = 10; % Minimum isolation distance
parameters.cluster.min_spikes    = 10;  % minimum number of spikes in cluster (change to mean firing rate threshold)

% Quality calculation parameters

parameters.cluster.thresh        = 4;   % number of standard deviations above background noise for mean amplitude of highest amplitude channel
parameters.cluster.thresh_xcorr  = 1;   % threshold for channel amplitude in order to be considered for cross-correlation


parameters.cluster.outlier_chi   = 0.0001; % threshold given as fraction of maximum of theoretical chi-squared function
parameters.cluster.outlier_std   = 3;      % max num. of SDs from centre of PCA cluster

parameters.cluster.pca_dims = 3; % first number of principle components to consider for spike filtering and cluster similarity calculation

% merge_thresh:
% Threshold to merge two clusters. Depends on 'parameters.sorting.method'
% - if method = KST, merge_thresh: sim score threshold (between 0 and 1)
parameters.cluster.merge_thresh = 0.85;

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