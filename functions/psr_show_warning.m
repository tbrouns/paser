function psr_show_warning(warnings)

nWarnings = length(warnings);
nLengths  = zeros(nWarnings,1);
for iWarning = 1:nWarnings
    nLengths(iWarning) = length(warnings{iWarning});
end

dashedline = repmat('-',1,max(nLengths));
col = '*blue';
cprintf(col,[dashedline '\n']); 
for iWarning = 1:nWarnings
    cprintf(col,[warnings{iWarning} '\n']);
end
cprintf(col,[dashedline '\n']);

end