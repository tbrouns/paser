function tf = psr_isfield(S,fname)

% PSR_ISFIELD - Check if field exists, allowing for non-existent parent fields
%
% Syntax:  tf = psr_isfield(S,fname)
%
% Inputs:
%    S     - Base structure
%    fname - String with full structure path to field we want to check
%
% Outputs:
%    tf - Returns true if field exists, false otherwise
%
% Example:
%    We want to know if the field "S.f1.f2" exists, so we
%    call: tf = psr_isfield(S,'S.f1.f2');

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

tf = false;
try 
    fname = split(fname,'.');
    fname = fname(2:end);
    for i = 1:length(fname)-1; S = [S.(fname{i})]; end
    fname = fname{end};
    tf = isfield(S,fname);
catch
    % Assume some parent-field doesn't exist
end

end