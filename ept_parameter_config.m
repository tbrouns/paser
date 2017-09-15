function parameters = ept_parameter_config()

% Hardware specifications
parameters.nelectrodes = 4; % Number of electrodes per array

% Spike detection

parameters.spikes.bp_high  = 6000;
parameters.spikes.bp_low   = 600;
parameters.spikes.bp_order = 10;

parameters.spikes.tsection = 60; % time of each individual section that is processed (in minutes)
parameters.cluster.method = 'fmm';

% Local field potential

parameters.lfp.bp_high  = 300;
parameters.lfp.bp_low   = 0.1;
parameters.lfp.bp_order =   5;

parameters.lfp.trial_length = 0.5; % sec

% Time frequency analysis

% see FT_FREQANALYSIS help

parameters.lfp.method = 'mtmconvol';
parameters.lfp.taper  = 'hanning';
parameters.lfp.pad    = 'nextpow2';

parameters.lfp.freq_lower = 2.0; % Hz
parameters.lfp.freq_upper = 60;  % Hz
parameters.lfp.freq_step  = 2.0; % Hz

parameters.lfp.time_step = 0.02; % sec

parameters.lfp.base_onset  = -0.3; % sec
parameters.lfp.base_offset = -0.1; % sec

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