function params = ept_analysis_parameters()

params.t_pre  = -50;
params.t_post = 100;
params.t_bin  = 5; % ms
params.t_array = [-400,0,50,100,200,400];

params.pad = 3; % 

params.smooth = false; 
params.sigma  =  25; % smoothing (ms)

end