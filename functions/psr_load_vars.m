function varargout = psr_load_vars(filepath,varnames)

% PSR_LOAD_VARS - Load variables from a MAT file, ignoring warnings
%
% Syntax:  varargout = psr_load_vars(filepath,varnames)
%
% Inputs:
%    filepath - Path to the MAT file
%    varnames - Variables we want to load, given as strings 
%
% Outputs:
%    varargout - Output variables
%
% Example:
%    [var_1,var_2] = psr_load_vars('myFile.mat',{'var_1','var_2'});

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

MSGID = 'MATLAB:load:variableNotFound';
warning('off', MSGID);
output = load(filepath,varnames{:});
nvars = length(varnames);
varargout = cell(nvars,1);
for iVar = 1:nvars
    varname = varnames{iVar};
    if (psr_isfield(output,['output.' varname]))
        varargout{iVar} = output.(varname);
    else
        varargout{iVar} = [];
    end
end
warning('on', MSGID);

end