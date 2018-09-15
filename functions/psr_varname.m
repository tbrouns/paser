function out = psr_varname(var)

% PSR_VARNAME - Returns name of input variable
%
% Syntax:  out = psr_varname(var)
%
% Inputs:
%    var - A variable
%
% Outputs:
%    out - Name of variable (string)
%
% Example:
%    var = 1;
%    out = psr_varname(var);
%    disp(out); % Prints 'var'

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

out = inputname(1);

end