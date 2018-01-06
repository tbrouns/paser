function params = psr_analysis_parameters()

params.Fs      = 30000;

%% Stimuli

params.stimOffset_1 = -37.1; % [ms]
params.stimOffset_2 = -15.1; % [ms]

%% LFP

params.lfp.window    = [-500,500];
params.lfp.base_win  = [-500,  0];
params.lfp.base_type = 'absolute';
params.lfp.base      = false;
params.lfp.plot_mean = false;

%% Spikes

params.t_bin   = 25; % [ms]
params.t_win   = [-600,600]; % Window to extract before and after stimulus [ms]
params.t_del   = [-0.5,0.5]; % window on both sides of stimulus to delete all spikes [ms]
params.t_pad   = 50; % [ms]
params.t_array = ...
    [-500, -50;  ...
    0,      50;  ...
    50,    100;  ...
    100,   200;  ...
    200,   400]; % [ms]

params.Nspikes = 10000;

% Baseline

params.base_win = [-500,-100];

% PSTH

params.psth_win = [-100,400]; % [ms]

% JPSTH

params.jpsth_win = [-100,400];

% Stimulus onset firing rate difference

params.diff_win = [50,100,150,200];

% ISI analysis

params.isi_bin = 1;   % ISI time bin [ms]
params.isi_max = 150; % Maximum ISI  [ms]

% ACF

params.acf_bin = 1;  % Bin size for spike binning [ms]
params.acf_max = 100; % [ms]
params.acf_win = 300; % [ms]

%

params.y_step_amp = 0.02; % max

params.smooth = false;
params.sigma  = 50; % smoothing [ms]

% Figures

params.fig_hght =  600;
params.fig_wdth = 1200;

end