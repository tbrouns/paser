function params = psr_analysis_parameters()

params.t_bin   = 5; % ms
params.t_win   = [-100,200]; % ms
params.t_array = [-400,0,50,100,200,400]; % ms
params.t_del   = 1; % window on both sides of stimulus to delete all spikes [ms]
params.pad = 3; % 

params.y_step_amp = 0.02; % max 

params.smooth = false; 
params.sigma  =  25; % smoothing (ms)

end