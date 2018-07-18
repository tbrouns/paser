function tf = isempty_field(S,fname)

% ISEMPTY_FIELD - Checks whether a field is empty/non-existent.
% This function is particularly useful when we want to check whether a
% deeply nested field is empty, but are not sure if the field and/or its
% parent fields exist. The function "isfield" cannot readily be applied for
% this purpose, because it requires the parent fields to exist.
% 
% Syntax:  tf = isempty_field(S,fname)
%
% Inputs:
%    S     - Structure or cell array of structures (vector)
%    fname - String with full structure path to field we want to check
%
% Outputs:
%    tf - Returns true if field is empty OR does not exist. Returns
%    false when field exists AND is not empty. Returns a logical vector
%    when input is cell array of structures.
% 
% Example:
%    We want to know if the field "S.f1.f2" is empty/non-existent, so we
%    call: tf = isempty_field(S,'S.f1.f2');
%
% Other m-files required: fieldnamesr

% Author: Terence Brouns 
% E-mail address: t.s.n.brouns@gmail.com 
% Date: 2018

if (iscell(fname)) % Input cell array of structures to check
    n  = length(fname); % Number of fields to check
    tf = true(n,1);
    for i = 1:n; tf(i) = isempty_field(S,fname{i}); end % Call function recursively
else
    tf    = true;
    fname = split(fname,'.');
    fname = fname(2:end);
    if (size(fname,2) > 0)
        fname = cell2mat(join(fname,'.'));
        if (~isempty(S))
            names = fieldnamesr(S,'full'); % Returns all fields
            I = strcmp(fname,names);
            if (any(I)) % field exists, now check if empty
                fname = split(fname,'.');
                for i = 1:length(fname); S = [S.(fname{i})]; end
                tf = isempty(S);
            end
        end
    end
end

end