function x = psr_gauss_smoothing(x,dt,sigma)

% PSR_GAUSS_SMOOTHING - Gaussian smoothing of time series data
%
% Syntax:  x = psr_gauss_smoothing(x,dt,sigma)
%
% Inputs:
%    x     - Time-series data
%    dt    - Sampling period
%    sigma - Standard deviation of Gaussian
%
% Outputs:
%    x - Gaussian smoothed time-series data
%
% See also: CONV,  GAUSSWIN

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

n = ceil(sigma / dt);
g = gausswin(n); % create Gaussian smoothing window
g = g / sum(g); % normalize
x = conv(x, g, 'same');

end