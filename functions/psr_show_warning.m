function psr_show_warning(warnings,isPATH,filename)

% PSR_SHOW_WARNING - Display a blue message in the command window
%
% Syntax:  psr_show_warning(warnings,isPATH)
%
% Inputs:
%    warnings - Cell array of strings containing messages we want to print
%               to command window
%    isPATH   - (Optional) Boolean indicating whether we want to display a path
%    filename - (Optional) String giving the filename where the warning
%               originated from. Typically you would want to use MFILENAME
%               as input here.
% 
% Example: 
%    psr_show_warning({'This','is','a','message'},false,mfilename)

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (nargin < 2); isPATH   = false; end
if (nargin < 3); filename = [];    end

if (~isempty(filename))
    warnings(2:end+1) = warnings;
    warnings{1} = ['ERROR in ' upper(filename) ':'];
end

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