function params = psr_analysis_parameters()

params.Fs = 30000;

params.t_bin   = 25; % [ms]
params.t_win   = [-1000,1000]; % Window to extract before and after stimulus [ms]
params.t_array = [-400,0,50,100,200,400]; % [ms]
params.t_del   = 1; % window on both sides of stimulus to delete all spikes [ms]
params.pad     = 3; % 

params.Nspikes = 10000;

% PSTH

params.psth_win = [-900,900]; % [ms]

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