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

%% General
parameters.general.precision = 1; % Number of digits after decimal point to keep for int16 data conversion

%% Development
parameters.develop.comparison = false;
parameters.develop.timing     = false;
parameters.develop.gamma      = 0.3;
parameters.develop.epsilon    = 0.5; % [ms]

%% PSR_ARTIFACT_FFT
% Parameters used in raw signal filtering in power spectrum

parameters.filter.fft.freq   = 10;   % Size of frequency window [Hz]
parameters.filter.fft.pad    = 0.01; % Size of half-window to extract artifact [Hz]
parameters.filter.fft.thresh = 5;    % Threshold to detect artifact peaks, given by number of STDs above mean (using MAD)

%% PSR_ARTIFACT_DIFF
parameters.filter.diff.win_slope    = 0.10; % [ms]
parameters.filter.diff.win_stimulus = 100;  % [ms]
parameters.filter.diff.win_pulse    = 30;   % Search window for second peak [ms]
parameters.filter.diff.win_artifact = 2.0;  % [ms]
parameters.filter.diff.win_padding  = 0.5;  % [ms]
parameters.filter.diff.thresh       = 10;  

%% SPIKE DETECTION

parameters.spikes.thresh      = 3.0; % Threshold given as number of standard deviations above mean for spike detection
parameters.spikes.ref_period  = 1.5; % Refractory period of spikes [ms]
parameters.spikes.window_size = 1.5; % Width of a spike [ms]. Take samples symmetrically around peak.
parameters.spikes.max_desync  = 0.25; % Maximum temporal desynchronization between channels of probe [ms]

% Zero-phase digital band-pass filter parameters for spike detection
parameters.spikes.bp.upper = 6000; % Upper cutoff frequency [Hz]
parameters.spikes.bp.lower = 600;  % Lower cutoff frequency [Hz]
parameters.spikes.bp.order = 10;   % Order of filter 

parameters.spikes.twin = 10; % Time of each individual section for spike detection (minutes)

%% SPIKE SORTING

% parameters.sorting.method
% 'CBP' : Continuous Basis Pursuit      - https://github.com/chinasaur/CBPSpikesortDemo
% 'FMM' : Focused Mixture Model         - https://github.com/decarlson/FMMSpikeSorter
% 'ISO' : ISO-SPLIT                     - https://github.com/magland/isosplit_old
% 'KFM' : Kalman filter mixture model   - https://github.com/AnaCalabrese/KFMM
% 'KST' : KiloSort                      - https://github.com/cortex-lab/KiloSort
% 'OPS' : OPASS                         - https://github.com/decarlson/opass
% 'OST' : OSort                         - http://www.urut.ch/new/serendipity/index.php?/pages/osort.html
% 'SPC' : Super-paramagnetic clustering - http://redishlab.neuroscience.umn.edu/MClust/MClust.html
% 'UMS' : UltraMegaSort2000             - https://github.com/danamics/UMS2K

parameters.sorting.method  = 'KST'; % Spike sorting method (see above)

%% Focused Mixture Model [FMM] (Dictionary learning) 
parameters.sorting.fmm.p     = 1e-4;  % Changes how aggresively to cluster, range 0-1 (0: less clustering, 1: more clustering)
parameters.sorting.fmm.k     = 5; 
parameters.sorting.fmm.kmax  = 16;    % Maximum number of clusters
parameters.sorting.fmm.align = false; % Whether to align the waveforms 

%% Continuous Basis Pursuit [CBP]
% [Need to set path to repository with "parameters.path.cbp"]
parameters.sorting.cbp.nC     = 16; % number of clusters

%% UltraMegaSort2000 [UMS]
parameters.sorting.ums.kmeans_size = 0.01; % target size for miniclusters as fraction of total number of spikes
parameters.sorting.ums.agg_cutoff  = 0.0001;  % higher = less aggregation, lower = more aggregation

%% Opass [OPS]
parameters.sorting.ops.dims    = 3;
parameters.sorting.ops.alph    = 1e-1;
parameters.sorting.ops.kappa_0 = 0.01;
parameters.sorting.ops.nu_0    = 0.1;
parameters.sorting.ops.Phi_0f  = 0.1; % Factor for 'Phi_0'
parameters.sorting.ops.a_pii   = 1;
parameters.sorting.ops.b_pii   = 1e7;
parameters.sorting.ops.betf    = 30; % Factor for 'bet'

%% Super paramagnetic clustering [SPC]
parameters.sorting.spc.mcs    = 0.01; % Minimum cluster size as fraction of total number of spikes
parameters.sorting.spc.dims   = 10;   % Number of inputs to the clustering 
parameters.sorting.spc.scales = 4;    % Number of scales for the wavelet decomposition 

%% ISO-SPLIT [ISO]
parameters.sorting.iso.mcs    = 0.01; % Minimum cluster size as fraction of total number of spikes
parameters.sorting.iso.dims   = 10;   % Number of inputs to the clustering 
parameters.sorting.iso.scales = 4;    % Number of scales for the wavelet decomposition 

%% Kalman filter mixture model [KFM]
parameters.sorting.kfm.nC     = 16; % number of clusters
parameters.sorting.kfm.dims   = 10; % Number of inputs to the clustering 
parameters.sorting.kfm.scales = 4;  % Number of scales for the wavelet decomposition 

%% KiloSort [KST]

% [Need to set path to repository with "parameters.path.kst"]
parameters.sorting.kst.GPU     = 1; % use CUDA (GPU computing)
parameters.sorting.kst.verbose = 1; % whether to print command line progress (true or false)

parameters.sorting.kst.Nfilt   = 32; % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)

% other options for controlling the model and optimization
parameters.sorting.kst.Nrank       = 3;        % matrix rank of spike template model (3)
parameters.sorting.kst.nfullpasses = 6;        % number of complete passes through data during optimization (6)
parameters.sorting.kst.maxFR       = 20000;    % maximum number of spikes to extract per batch (20000)
parameters.sorting.kst.ntbuff      = 64;       % samples of symmetrical buffer for whitening and spike detection
parameters.sorting.kst.NT          = 128*1024+ parameters.sorting.kst.ntbuff; % this is the batch size (try decreasing if out of memory)
% for GPU should be multiple of 32 + ntbuff

% the following options can improve/deteriorate results.
% when multiple values are provided for an option, the first two are beginning and ending anneal values,
% the third is the value used in the final pass.
parameters.sorting.kst.Th               = [6 12 12];    % threshold for detecting spikes on template-filtered data ([6 12 12])
parameters.sorting.kst.lam              = [10 30 30];   % large means amplitudes are forced around the mean ([10 30 30])
parameters.sorting.kst.nannealpasses    = 4;            % should be less than nfullpasses (4)
parameters.sorting.kst.momentum         = 1./[20 1000]; % start with high momentum and anneal (1./[20 1000])
parameters.sorting.kst.shuffle_clusters = 1;            % allow merges and splits during optimization (1)
parameters.sorting.kst.mergeT           = .1;           % upper threshold for merging (.1)
parameters.sorting.kst.splitT           = .1;           % lower threshold for splitting (.1)

% options for initializing spikes from data
parameters.sorting.kst.initialize      = 'fromData';    %'fromData' or 'no'
parameters.sorting.kst.spkTh           = -4;      % spike threshold in standard deviations (4)
parameters.sorting.kst.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
parameters.sorting.kst.long_range      = [30 6];  % ranges to detect isolated peaks ([30 6])
parameters.sorting.kst.maskMaxChannels = 5;       % how many channels to mask up/down ([5])
parameters.sorting.kst.crit            = .65;     % upper criterion for discarding spike repeates (0.65)
parameters.sorting.kst.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)

%% CLUSTER PARAMETERS & QUALITY CONTROL

%% PSR_SST_CLUSTER_CLASSIFIER

% Hard thresholds

parameters.cluster.max_rpv       = 0.10; % Maximum fraction of refractory period violations (RPVs) in cluster
parameters.cluster.max_sub       = 0.10; % Maximum fraction of sub-threshold spikes in cluster
parameters.cluster.max_amplitude = 200;  % Maximum absolute mean spike amplitude of cluster [microvolt]
parameters.cluster.max_p2p       = 300;  % Maximum absolute mean peak-to-peak amplitude of cluster [microvolt]
parameters.cluster.min_spikes    = 50;   % Minimum number of spikes in cluster
parameters.cluster.min_frate     = 0.05; % Minimum average firing rate (Hz) 
parameters.cluster.max_lratio    = 0.25; % Maximum L-ratio for single unit
parameters.cluster.min_isodist   = 20;   % Minimum isolation distance for single unit

%% PSR_SST_FILTER_SPIKES

parameters.filter.spikes.corr_global = true; % Remove high global correlation spikes
parameters.filter.spikes.mse_cluster = true; % Remove spikes with low within-cluster correlations
parameters.filter.spikes.rpvs        = true; % Remove refractory period violations
parameters.filter.spikes.amp         = true; % Remove spikes with amplitude at wrong channel or position

parameters.filter.spikes.corr_thresh_global = 0.5; % Correlation across probes threshold for spike removal
parameters.filter.spikes.corr_thresh_chan   = 0.5; % Used for non-spike channel noise subtraction
parameters.filter.spikes.mse_thresh         = 5.0; % 
parameters.filter.spikes.amp_offset         = 0.2; % Normalized amplitude difference from maximum, for channel selection

% Quality calculation parameters

parameters.cluster.amplitude_nbin = 20; % Average number of spikes per bin 

parameters.cluster.thresh       = 4; % Number of standard deviations above mean for average amplitude of highest amplitude channel
parameters.cluster.thresh_xcorr = 1; % Threshold for channel amplitude in order to be considered for cross-correlation, given by number of STDs.

%% PSR_SST_AMP_SPLIT & 

parameters.cluster.split.bin      = 20;   % Average number of spikes per bin
parameters.cluster.split.smooth   = 7;    % Number of bins for Gaussian smoothing
parameters.cluster.split.thresh   = 100;  % Number of spikes
parameters.cluster.split.prom     = 0.1;  % Minimum fraction of maximum peak
parameters.cluster.split.width    = 3;    % Minimum width of neighbouring peaks 
parameters.cluster.split.terminal = 0.05; % Threshold for termination (fraction of maximum peak)

% Cluster merge criterion
parameters.cluster.merge_thresh = 0.90; % Cross-correlation threshold of concatenated channels

% Stability 
parameters.cluster.stability_win = 10; % [sec]

% Mixture of drifting t-distribution model [MDT]
parameters.cluster.mdt.dims              = 10; % Number of dimensions for wavelet decomposition                   
parameters.cluster.mdt.scales            = 4;  % Number of scales for the wavelet decomposition 
parameters.cluster.mdt.nu                = 7;  % t-distribution nu parameter (smaller = heavier tails) 
parameters.cluster.mdt.q_perhour         = 2;  % Drift regularization (smaller = more smoothing)         
parameters.cluster.mdt.timeframe_minutes = 1;  % Time frame duration (mostly a computational thing)   

%% LOCAL FIELD POTENTIAL

parameters.lfp.rsfactor = 4; % Multiple of high band-pass filter frequency to resample at

parameters.lfp.bp.upper = 300; % Band-pass filter upper cut-off frequency [Hz]
parameters.lfp.bp.lower = 0.1; % Band-pass filter lower cut-off frequency [Hz]
parameters.lfp.bp.order =   5; % Band-pass filter order

parameters.lfp.trial.onset  = -0.2; % sec
parameters.lfp.trial.offset =  0.5; % sec
parameters.lfp.trial.padding = 1.0; % time window before and after trial on- and offset (sec)

parameters.lfp.miss_thresh = 0.5; % Maximum fraction of missing data in fft window

%% FIELDTRIP PARAMETERS
% [Need to set path to repository with "parameters.path.ft"]

% For details on the following parameters, see FT_FREQANALYSIS help.

parameters.lfp.method = 'mtmconvol';
parameters.lfp.taper  = 'hanning';
parameters.lfp.pad    = 'nextpow2';

parameters.lfp.freq.lower = 2.0; % Lower limit of frequency bins [Hz]
parameters.lfp.freq.upper = 60;  % Upper limit of frequency bins [Hz]
parameters.lfp.freq.step  = 0.5; % Frequency bin size [Hz]

parameters.lfp.time_step = 0.02; % Time bin size [sec]

parameters.lfp.ncycles = 5; % number of cycles for frequency in sliding time window

%% Local field potential - artifact removal 

% Power spectral density (PSD)
parameters.lfp.artifact.freqRange  = linspace(10,100,64); % Frequencies to calculate PSD for
parameters.lfp.artifact.threshPSD  = 3; % PSD threshold
parameters.lfp.artifact.window     = 0.5; % Window to calculate PSD [sec]

% Amplitude
parameters.lfp.artifact.cat            = true;  % Concatenate sessions
parameters.lfp.artifact.sigma          = 0.005; % Gaussian filter window [sec]
parameters.lfp.artifact.tsection       = 10;    % Window length in which we detect background noise [sec]
parameters.lfp.artifact.threshAmpUpper =  6;    % Number of standard deviations above background noise
parameters.lfp.artifact.threshAmpLower =  2; 
parameters.lfp.artifact.tSlope         =  6;    % [ms]

% General
parameters.lfp.artifact.interval = 0.50; % Minimum clean interval between artifacts [sec]
parameters.lfp.artifact.padding  = 0;    % Padding added on either side of detected artifact [sec]

%% Stimulus onset

%% PSR_MFA_DETECTION
% Magnetic field artifact (MFA) parameters

parameters.mfa.threshold = 0.99; % Absolute signal threshold
parameters.mfa.min_dur   = 0.1; % Minimum duration of pulse [sec]

%% PSR_CAM_DETECTION

parameters.cam.thresh = 10;

%% Experiment specific

% Active vs. passive touch
parameters.experimental.ap.process = false; % Process A.vs.P. sections
parameters.experimental.ap.diff    = true;
parameters.experimental.ap.offset  = true;