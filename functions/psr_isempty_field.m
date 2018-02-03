function tf = psr_isempty_field(S,fname)

% Returns true if field is empty OR doesn't exist. Returns false when
% field exists AND is not empty.
% 
% X: structure
% fname: string of field of interest
% 
% Example: 
%
% We want to know if the field "S.F1.F2" is empty or not. We call the
% function as follows: I = psr_isempty_field(S,'F1.F2');

tf = true;
if (~isempty(S))
    names = [fieldnames(S);fieldnamesr(S)];
    I = strcmp(fname,names);
    if (any(I)) % field exists, now check if empty
        fname = split(fname,'.');
        for i = 1:length(fname); S = S.(fname{i}); end
        tf = isempty(S);
    end
end

end