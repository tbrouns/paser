function data = psr_lfp_artifact_removal(data,artifacts,parameters)    

% Signal can be given as single time-series (Nchans x Npoints), or as a
% cell array of such matrices or FieldTrip data structures

if (iscell(data)); nTrials = length(data);
else,              nTrials = 1;
end

% Set constants

fields   = fieldnames(artifacts);
nFields  = length(fields);
Fs       = parameters.Fr; 
interval = Fs * parameters.lfp.artifact.interval;
padding  = Fs * parameters.lfp.artifact.padding;

artifactsAll = [];

for iTrial = 1:nTrials
    
    dataTrial = data{iTrial};
    if (isfield(dataTrial,'trial')); dataTrial = dataTrial.trial{1}; end
    slength = size(dataTrial,2);
    
    for iField = 1:nFields
        
        artifactsTrial = artifacts.(fields{iField}){iTrial};
        
        if (~isempty(artifactsTrial))
        
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

        end
        
        artifactsAll.(fields{iField}){iTrial} = artifactsTrial;

        
    end
    
    if (isfield(data{iTrial},'trial')); data{iTrial}.trial = {dataTrial};
    else,                               data{iTrial}       =  dataTrial;
    end
end

% Data conversion

for iTrial = 1:nTrials
    if (isfield(data{iTrial},'trial'))
        for iField = 1:nFields
            data{iTrial}.artifacts.(fields{iField}) = artifactsAll.(fields{iField}){iTrial};
        end
    end
end

end