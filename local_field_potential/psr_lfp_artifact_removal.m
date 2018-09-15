function data = psr_lfp_artifact_removal(data,artifacts,parameters)

% PSR_LFP_ARTIFACT_REMOVAL - Remove detected artifacts from data 
% Set data to NaN at locations of the detected input artifacts. 
%  
% Syntax:  data = psr_lfp_artifact_removal(data,artifacts,parameters)
%
% Inputs:
%    data       - Cell array of time-series data. Each cell contains a
%                 FieldTrip data structure or an array of shape:
%                 [Number of channels x Number of data points]
%    artifacts  - Two-column array, where the first column are the
%                 onsets and the second column the offsets of the detected
%                 power spectral density artifacts. Timestamps given in
%                 secs.
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    data - Time-series data where artifact regions have been set to NaN
%
% See also: PSR_LFP_ARTIFACT_DETECTION_AMP, PSR_LFP_ARTIFACT_DETECTION_PSD

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Set constants
nBlocks  = length(data);
fields   = fieldnames(artifacts);
nFields  = length(fields);
Fs       = parameters.Fr; 
interval = Fs * parameters.lfp.artifact.interval;
padding  = Fs * parameters.lfp.artifact.padding;

artifactsAll = [];

for iBlock = 1:nBlocks
    
    dataBlock = data{iBlock};
    if (isfield(dataBlock,'trial')); dataBlock = dataBlock.trial{1}; end
    sLength = size(dataBlock,2);
    
    for iField = 1:nFields
        
        artifactsBlock = artifacts.(fields{iField}){iBlock};
        
        if (~isempty(artifactsBlock))
        
            artifactsBlock = round(Fs * artifactsBlock) + 1;
            artifactsBlock = sortrows(artifactsBlock);
            onsets  = artifactsBlock(:,1);
            offsets = artifactsBlock(:,2);

            % Add padding
            onsets  = onsets  - padding;
            offsets = offsets + padding;

            % Combine artifact intervals
            d  = onsets(2:end) - offsets(1:end-1);
            del = find(d <= interval); % adjacent intervals
            onsets (del + 1) = [];
            offsets(del)     = [];

            % Check limits
            onsets ( onsets < 1)       = 1;
            offsets(offsets > sLength) = sLength;
            artifactsBlock = [onsets,offsets];
            artifactsBlock = unique(artifactsBlock,'rows');

            nArtifacts = size(artifactsBlock,1);
            for iArtifact = 1:nArtifacts
                x = artifactsBlock(iArtifact,1):artifactsBlock(iArtifact,2);
                dataBlock(:,x) = NaN;
            end
        end
        
        artifactsAll.(fields{iField}){iBlock} = artifactsBlock;
    end
    
    if (isfield(data{iBlock},'trial')); data{iBlock}.trial = {dataBlock};
    else,                               data{iBlock}       =  dataBlock;
    end
end

% Data conversion

for iBlock = 1:nBlocks
    if (isfield(data{iBlock},'trial'))
        for iField = 1:nFields
            data{iBlock}.artifacts.(fields{iField}) = artifactsAll.(fields{iField}){iBlock};
        end
    end
end

end