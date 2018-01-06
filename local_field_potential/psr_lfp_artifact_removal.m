function [data,artifacts] = psr_lfp_artifact_removal(data,artifacts,parameters)    

% Set constants
nTrials  = size(artifacts,1);
nTypes   = size(artifacts,2);
Fs       = parameters.Fr; 
interval = Fs * parameters.lfp.artifact.interval;
padding  = Fs * parameters.lfp.artifact.padding;

artifactsAll = cell(nTrials,1);

for iTrial = 1:nTrials
    
    dataTrial = data{iTrial};
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
    onsets(  onsets < 1)       = 1;
    offsets(offsets > slength) = slength;
    artifactsTrial = [onsets,offsets];
    artifactsTrial = unique(artifactsTrial,'rows');
        
    nArtifacts = size(artifactsTrial,1);
    for iArtifact = 1:nArtifacts
        x = artifactsTrial(iArtifact,1):artifactsTrial(iArtifact,2);
        dataTrial(:,x) = NaN;
    end
    
    artifactsAll{iTrial} = artifactsTrial;
    data{iTrial}         = dataTrial;
end

artifacts = artifactsAll; % Output

end