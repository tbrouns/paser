function tf = psr_isfield(S,fname)

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