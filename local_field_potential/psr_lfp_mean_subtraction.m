function freq = psr_lfp_mean_subtraction(freq)

% PSR_LFP_MEAN_SUBTRACTION - Common average referencing
%
% Syntax:  freq = psr_lfp_mean_subtraction(freq)
%
% Inputs:
%    freq - FieldTrip LFP data structure (see README)
%
% Outputs:
%    freq - Data after subtracting mean across all channels
%
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Subtract channel mean from all probe channels
% Common average referencing 

nBlocks = length(freq);
for iBlock = 1:nBlocks
    dataBlock = freq{iBlock};
    if (isfield(dataBlock,'trial'))
        dataBlock = dataBlock.trial{1};
        nChans = size(dataBlock,1);
        if (nChans > 1); freq{iBlock}.trial{1} = meanSubtraction(dataBlock); end
    end
end

end

function data = meanSubtraction(data)

dataMean = nanmean(data,1); % Mean across channels
data = bsxfun(@minus,data,dataMean);

end