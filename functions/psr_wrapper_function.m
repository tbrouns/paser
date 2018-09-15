function output = psr_wrapper_function(cfg)

% PSR_WRAPPER_FUNCTION - Function for calling custom analysis function
%
% Syntax:  output = psr_wrapper_function(cfg)
%
% Inputs:
%    cfg - Structure that should at least contain an "fpath" field that
%          points to the custom analysis function. The "cfg" structure is
%          also given as input to the custom analysis function.
%
% Outputs:
%    output - Depends on the custom analysis function
%
% See also: PSR_BATCH_ANALYSIS

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

fpath = cfg.fpath; % Path to custom analysis function
names = split(fpath,'\');
names = names(~cellfun('isempty',names)); % remove empty cells
fname = names{end};
fpath = join(names(1:end-1),'\');

oldFolder = cd(fpath{1}); % cd to path to ensure we call correct function
output = feval(fname,cfg); % Call the custom analysis function
cd(oldFolder); % go back to where we came from

end