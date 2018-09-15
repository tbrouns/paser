function x = psr_int16_to_single(x,parameters)

% PSR_INT16_TO_SINGLE - Convert int16 array to single precision
%
% Syntax:  x = psr_int16_to_single(x,parameters)
%
% Inputs:
%    x          - 16-bit signed integer array
%    parameters - See README
%
% Outputs:
%    x - Single precision array

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (isa(x,'int16'))
    precision = 10^parameters.general.precision;
    x = single(x) / precision;
end
end