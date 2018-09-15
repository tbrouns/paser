function [data,nBlocks] = psr_lfp_conversion(data)

% PSR_LFP_CONVERSION - Convert local field potential data format
%
% Syntax:  [data,nBlocks] = psr_lfp_conversion(data)
%
% Inputs:
%    data - Can be given as: 
%           (1) Single time-series with shape: [Number of channels x Number of data points]
%           (2) A cell array of the matrix type given by (1) 
%           (3) Cell array of FieldTrip data structures
%
% Outputs:
%    data    - Data in format given by (2)
%    nBlocks - Number of experimental blocks

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (iscell(data))
    nBlocks = length(data);
else % Single time-series
    data = {data};
    nBlocks = 1;
end

% Check for FieldTrip, and do conversion if necessary

for iBlock = 1:nBlocks
    if isfield(data{iBlock},'trial') 
        data{iBlock} = data{iBlock}.trial{1};
    end
end

end

%------------- END OF CODE --------------