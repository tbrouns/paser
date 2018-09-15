function [y,x] = psr_remove_field(y,str)

% PSR_REMOVE_FIELD - Removes field from structure and returns removed field
%
% Syntax:  [y,x] = psr_remove_field(y,str)
%
% Inputs:
%    y   - Structure for which we want to remove a field
%    str - Name of field we want to remove (string)
%
% Outputs:
%    y - Structure with removed field
%    x - Data from removed field

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Removes field and returns removed field
try
    x = y.(str);
    y = rmfield(y,str);
catch 
    % Assumes that field does not exist
    x = [];
end

end