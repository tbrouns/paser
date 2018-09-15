function k = psr_session_strcmp(parameters,str)

% PSR_SESSION_STRCMP - Check which sessions match the input string
%
% Syntax:  k = psr_session_strcmp(parameters,str)
%
% Inputs:
%    parameters - See README
%    str        - Input string 
%
% Outputs:
%    k - Logical array with length equal to the number of sessions.
%        Elements corresponding to sessions that contain the input string
%        are set to true

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

nSessions = size(parameters.general.session,2);
sessions  = cell(nSessions,1);
for iSession = 1:nSessions % Convert to lower case
    sessions{iSession} = lower(parameters.general.session{iSession});
end

% Find active and passive trials
I = find(~cellfun(@isempty,strfind(sessions,str)));
k = false(size(parameters.general.sessionIndex));
k(ismember(parameters.general.sessionIndex,I)) = true;

end