function [data,artifacts] = psr_lfp_artifact_removal(data,artifacts,parameters)    

% Signal can be given as single time-series (Nchans x Nsamples), or as a
% cell array of such matrices or FieldTrip data structures

if (iscell(data)); nTrials = length(data);
else,              nTrials = 1;
end

% Set constants

nTypes   = size(artifacts,2);
Fs       = parameters.Fr; 
interval = Fs * parameters.lfp.artifact.interval;
padding  = Fs * parameters.lfp.artifact.padding;

artifactsAll = cell(nTrials,1);

for iTrial = 1:nTrials
    
    dataTrial = data{iTrial};
    if (isfield(dataTrial,'trial')); dataTrial = dataTrial.trial{1}; end
    slength = size(dataTrial,2);
    artifactsTrial = [];
    for iType = 1:nTypes; artifactsTrial = [artifactsTrial;artifacts{iTrial,iType}]; end %#ok
    artifactsTrial = round(Fs * artifactsTrial) + 1;
    artifactsTrial = sortrows(artifactsTrial);
    onsets  = artifactsTrial(:,1);
    offsets = artifactsTrial(:,2);
    
    % Add padding
    onsets  = onsets  - padding;
    offsets = offsets + padding;
    
    % Combine artifact intervals
    d  = onsets(2:end) - offsets(1:end-1);
    id = find(d <= interval); % adjacent intervals
    onsets (id + 1) = [];
    offsets(id)     = [];
    
    % Check limits
    onsets ( onsets < 1)       = 1;
    offsets(offsets > slength) = slength;
    artifactsTrial = [onsets,offsets];
    artifactsTrial = unique(artifactsTrial,'rows');
    
    nArtifacts = size(artifactsTrial,1);
    for iArtifact = 1:nArtifacts
        x = artifactsTrial(iArtifact,1):artifactsTrial(iArtifact,2);
        dataTrial(:,x) = NaN;
    end
        
    artifactsAll{iTrial} = artifactsTrial;
    
    if (isfield(data{iTrial},'trial')); data{iTrial}.trial = {dataTrial};
    else,                               data{iTrial}       =  dataTrial;
    end
    
end

for iTrial = 1:nTrials
    if (isfield(data{iTrial},'trial'))
        data{iTrial}.artifacts = artifactsAll{iTrial};
    end
end

artifacts = artifactsAll; % Output

end