% PSR_PARAMETER_DEFAULT - Script to set default parameters in PASER.
% See comments before and after each parameter for its purpose.
%
% Syntax:  psr_parameter_default
%
% Outputs:
%    parameters - Structure that contains parameters for data processing
%

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

%% Pipeline conditional parameters
% Used to skip certain sections of data processing pipeline

parameters.process.spikes = true; % Perform spike detection
parameters.process.lfp    = true; % Perform LFP detection
parameters.process.active_passive = true; % Active vs. passive conditions

% General

parameters.general.savelist = {'spikes','metadata','freq','parameters'}; % What variables to save in output MAT file
parameters.general.precision = 1; % decimal place to round raw int16 data to (should be integer)

%% PSR_ARTIFACT_FFT
% Parameters used in raw signal filtering in power spectrum

parameters.filter.fft_freq   = 10;   % Size of frequency window [Hz]
parameters.filter.fft_pad    = 0.01; % Size of half-window to extract artifact [Hz]
parameters.filter.fft_thresh = 5;    % Threshold to detect artifact peaks, given by number of STDs above mean (using MAD)

%% SPIKE DETECTION

% parameters.spikes.method:
% 'auto'    : threshold calculated from background noise
% 'manual'  : user defined threshold 
% 'mad'     : median absolute deviation [MAD] (default)

parameters.spikes.method = 'mad'; % Spike detection method (see above)

% parameters.spikes.thresh:
% If 'auto'   : number of standard deviations above mean
% If 'manual' : absolute threshold in microvolts
% If 'mad'    : number of standard deviations above mean using median absolute deviation

parameters.spikes.thresh      = 3.0; % Threshold for spike detection (see above)

parameters.spikes.ref_period  = 1.5; % Refractory period of spikes [ms]
parameters.spikes.window_size = 1.5; % Width of a spike [ms]. Take samples symmetrically around peak.
parameters.spikes.max_desync  = 0.2; % Maximum temporal desynchronization between channels of probe [ms]

parameters.spikes.artifacts_corr = 0.5; % Minimum correlation threshold for spike artifact detection across tetrodes

% Zero-phase digital band-pass filter parameters for spike detection
parameters.spikes.bp_upper = 6000; % Upper cutoff frequency [Hz]
parameters.spikes.bp_lower = 600;  % Lower cutoff frequency [Hz]
parameters.spikes.bp_order = 10;   % Order of filter 

parameters.spikes.twin = 10; % Time of each individual section for spike detection (minutes)

parameters.spikes.artifacts_removal = true;  % Remove artifacts based on correlation across tetrodes

%% SPIKE SORTING

% parameters.sorting.method
% 'CBP' : Continuous Basis Pursuit      - https://github.com/chinasaur/CBPSpikesortDemo
% 'FMM' : Focused Mixture Model         - https://github.com/decarlson/FMMSpikeSorter
% 'ISO' : ISO-SPLIT                     - https://github.com/magland/isosplit_old
% 'KST' : KiloSort                      - https://github.com/cortex-lab/KiloSort
% 'OPS' : OPASS                         - https://github.com/decarlson/opass
% 'SPC' : Super-paramagnetic clustering - http://redishlab.neuroscience.umn.edu/MClust/MClust.html
% 'UMS' : UltraMegaSort2000             - https://github.com/danamics/UMS2K

parameters.sorting.method  = 'KST'; % Spike sorting method (see above)

% Focused Mixture Model [FMM] (Dictionary learning) 
parameters.sorting.fmm.p     = 1e-4;  % Changes how aggresively to cluster, range 0-1 (0: less clustering, 1: more clustering)
parameters.sorting.fmm.k     = 5; 
parameters.sorting.fmm.kmax  = 16;    % Maximum number of clusters
parameters.sorting.fmm.align = false; % Whether to align the waveforms 

% Continuous Basis Pursuit [CBP]
% [Need to set path to repository with "parameters.path.cbp"]

% UltraMegaSort2000 [UMS]
parameters.sorting.ums.kmeans_size = 0.01; % target size for miniclusters as fraction of total number of spikes
parameters.sorting.ums.agg_cutoff  = 0.0001;  % higher = less aggregation, lower = more aggregation

% Super paramagnetic clustering [SPC]
parameters.sorting.spc.dims = 3;

% KiloSort [KST]
% [Need to set path to repository with "parameters.path.kst"]
parameters.sorting.kst.cuda = 1; % use CUDA (GPU computing)

%% CLUSTER PARAMETERS & QUALITY CONTROL

% Thresholds

parameters.cluster.max_rpv       = 0.05; % Maximum fraction of refractory period violations (RPVs) in cluster
parameters.cluster.max_missing   = 0.05; % Maximum fraction of missing spikes in cluster due to threshold
parameters.cluster.max_amplitude = 300;  % Maximum absolute mean spike amplitude of cluster (microvolt)
parameters.cluster.max_artifact  = 1.25; % Maximum ratio between actual and expected number of spikes in LFP artifact region
parameters.cluster.max_lratio    = 1.0;  % Maximum L-ratio for single unit
parameters.cluster.max_lag       = parameters.spikes.max_desync; % Maximum lag of maximum pairwise cross correlation (ms)

parameters.cluster.min_isodist   = 10; % Minimum isolation distance for single unit
parameters.cluster.min_spikes    = 10; % Minimum number of spikes in cluster (change to mean firing rate threshold)

% Quality calculation parameters

parameters.cluster.thresh       = 4; % Number of standard deviations above mean for average amplitude of highest amplitude channel
parameters.cluster.thresh_xcorr = 1; % Threshold for channel amplitude in order to be considered for cross-correlation, given by number of STDs.

parameters.cluster.outlier_chi = 0.0001; % Threshold given as fraction of maximum of theoretical chi-squared function
parameters.cluster.outlier_std = 3;      % Maximum number of STDs from centre of PCA cluster

parameters.cluster.pca_dims = 3; % first number of principle components to consider for spike filtering and cluster similarity calculation

% parameters.cluster.merge_thresh:
% Threshold to merge two clusters. Depends on 'parameters.sorting.method'
% If parameters.sorting.method = KST, then merge_thresh: sim score threshold (between 0 and 1)
parameters.cluster.merge_thresh = 0.85;

% Stability 

parameters.cluster.stability_win = 10; % [sec]

% Mixture of drifting t-distribution model [MDT]
parameters.cluster.mdt.dims              = 12; % Number of PC dimensions                               (default: 12)
parameters.cluster.mdt.nu                = 7;  % t-distribution nu parameter (smaller = heavier tails) (default: 7)
parameters.cluster.mdt.q_perhour         = 2;  % Drift regularization (smaller = more smoothing)       (default: 2)  
parameters.cluster.mdt.timeframe_minutes = 1;  % Time frame duration (mostly a computational thing)    (default: 1)

%% LOCAL FIELD POTENTIAL

parameters.lfp.rsfactor = 4; % Multiple of high band-pass filter frequency to resample at

parameters.lfp.bp_upper = 300;
parameters.lfp.bp_lower = 0.1;
parameters.lfp.bp_order =   5;

parameters.lfp.trial_onset  = -0.2; % sec
parameters.lfp.trial_offset =  0.5; % sec
parameters.lfp.trial_padding = 1.0; % time window before and after trial on- and offset (sec)

%% FIELDTRIP PARAMETERS
% [Need to set path to repository with "parameters.path.ft"]

% For details on the following parameters, see FT_FREQANALYSIS help.

parameters.lfp.method = 'mtmconvol';
parameters.lfp.taper  = 'hanning';
parameters.lfp.pad    = 'nextpow2';

parameters.lfp.freq_lower = 2.0; % Lower limit of frequency bins [Hz]
parameters.lfp.freq_upper = 60;  % Upper limit of frequency bins [Hz]
parameters.lfp.freq_step  = 0.5; % Frequency bin size [Hz]

parameters.lfp.time_step = 0.02; % Time bin size [sec]

parameters.lfp.ncycles = 5; % number of cycles for frequency in sliding time window

parameters.lfp.base_onset  = -0.3; % sec
parameters.lfp.base_offset = -0.1; % sec

parameters.lfp.base_type = 'relative';

parameters.lfp.artifact_freq     = 30; % Maximum artifact frequency [Hz]
parameters.lfp.artifact_tsection = 10; % Window length in which we detect background noise [sec]
parameters.lfp.artifact_thresh   =  7; % Number of standard deviations above background noise
parameters.lfp.artifact_frac     = .1; % Fraction of STD for determining artifact on- and offset based on derivative 
parameters.lfp.artifact_padding  =  1; % Padding added on either side of detected artifact [sec]

%% Stimulus onset

%% PSR_MFA_DETECTION
% Magnetic field artifact (MFA) parameters

parameters.mfa.p2ptime  = 0.25; % Time between min and max peaks (ms)
parameters.mfa.bp_upper = 6000; % Upper cutoff frequency [Hz]
parameters.mfa.bp_lower = 1000; % Lower cutoff frequency [Hz]
parameters.mfa.bp_order = 10;   % Filter order
parameters.mfa.thresh   = 8;    % Number of STDs above the mean using median absolute deviation
parameters.mfa.control  = 100;  % Number of random stimulus onset times for control data
parameters.mfa_combine = false; % Use both ADC and CONTINUOUS files to detect artifacts (unreliable)

%% PSR_CAM_DETECTION

parameters.cam.thresh = 10;