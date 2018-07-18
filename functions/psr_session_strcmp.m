function k = psr_session_strcmp(parameters,str)

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