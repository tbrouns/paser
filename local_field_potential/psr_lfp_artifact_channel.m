function data = psr_lfp_artifact_channel(data,parameters)

% PSR_LFP_ARTIFACT_CHANNEL - Remove high noise channels
%
% Syntax:  data = psr_lfp_artifact_channel(data,parameters)
%
% Inputs:
%    data       - Same as input for PSR_LFP_CONVERSION
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    data - Data without high noise channels
% 
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

dataAbs = psr_lfp_conversion(data);
dataAbs = cat(2,dataAbs{:}); % Combine data across all blocks

nChans  =    size(dataAbs,1);
dataAbs =     abs(dataAbs);
dataAvg = nanmean(dataAbs,2); 
dataStd =  nanstd(dataAbs,[],2);

chanDiff = zeros(nChans,1);

for iChan = 1:nChans
   
    I = true(nChans,1);
    I(iChan) = false;
    
    dataStdChan = nanmean(dataStd(I));
    dataAvgMean = nanmean(dataAvg(I));
    
    chanDiff(iChan) = (dataAvg(iChan) - dataAvgMean) / dataStdChan;
    
end

% Only keep channels that aren't all NaNs
keep = ~all(isnan(dataAbs),2);
data = removeChannels(data,keep);
chanDiff = chanDiff(keep);

% Only keep channels that have low channel difference
if length(chanDiff) > 1
    keep = chanDiff <= parameters.lfp.artifact.chan.thresh;
    data = removeChannels(data,keep);
end

end

function data = removeChannels(data,keep)
    
nBlocks = length(data);
for iBlock = 1:nBlocks
    dataBlock = data{iBlock};
    if (isfield(dataBlock,'trial'))
        data{iBlock}.trial{1} = dataBlock.trial{1}(keep,:);
        data{iBlock}.label    = dataBlock.label(keep);
    end
end

end

%------------- END OF CODE --------------