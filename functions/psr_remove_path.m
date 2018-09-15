function psr_remove_path(str)

% PSR_REMOVE_PATH - Removes paths containing an input string from search path
%
% Syntax:  psr_remove_path(str)
%
% Inputs:
%    str - Input string that is contained in the paths that we want to
%          remove from the search path
%
% See also: PSR_FT_PATH

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

paths_all = path;
paths_all = strsplit(paths_all,';');
k = ~cellfun(@isempty,strfind(paths_all,str)); %#ok
paths_del = paths_all(1,k);
paths_del = strjoin(paths_del,';');
rmpath(paths_del);

end