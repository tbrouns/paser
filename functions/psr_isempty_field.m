function tf = psr_isempty_field(S,fname)

% Returns true if field is empty OR doesn't exist. Returns false when
% field exists AND is not empty.
%
% X: structure
% fname: string of field of interest or cell array of strings
%
% Example:
%
% We want to know if the field "S.F1.F2" is empty or not. We call the
% function as follows: I = psr_isempty_field(S,'S.F1.F2');

if (iscell(fname))
    n  = length(fname); % number of fields to check
    tf = true(n,1);
    for i = 1:n; tf(i) = psr_isempty_field(S,fname{i}); end
else
    tf    = true;
    fname = split(fname,'.');
    fname = fname(2:end);
    if (size(fname,2) > 0)
        fname = join(fname,'.');
        fname = fname{1};
        if (~isempty(S))
            try
                names = fieldnamesr(S,'full');
                I = strcmp(fname,names);
                if (any(I)) % field exists, now check if empty
                    fname = split(fname,'.');
                    for i = 1:length(fname); S = S.(fname{i}); end
                    tf = isempty(S);
                end
            catch
            end
        end
    end
end

end