function psr_show_warning(warnings,isPATH)

if (nargin < 2); isPATH = false; end

nWarnings = length(warnings);

% Change backslash to forward slash to avoid errors
if (isPATH)
    for iWarning = 1:nWarnings
        str = warnings{iWarning};
        k = strfind(str,'\'); str(k) = '/';
        warnings{iWarning} = str;
    end
end

% Find longest string

nLengths  = zeros(nWarnings,1);
for iWarning = 1:nWarnings
    nLengths(iWarning) = length(warnings{iWarning});
end

% Print the error
dashedline = repmat('-',1,max(nLengths));
col = '*blue';
cprintf(col,[dashedline '\n']); 
for iWarning = 1:nWarnings
    cprintf(col,[warnings{iWarning} '\n']);
end
cprintf(col,[dashedline '\n']);

end