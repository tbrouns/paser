function combis = psr_ft_combis(spikesFT,unitIDs)

% PSR_FT_COMBIS - Return all possible unit combinations
%
% Syntax:  combis = psr_ft_combis(spikesFT,unitIDs)
%
% Inputs:
%    spikesFT - A FieldTrip spike structure
%    unitIDs  - Unit IDs we want to find all the combinations for
% 
% Outputs:
%    combis - List of all possible unit combinations, with shape:
%             [Number of combinations x 2]
%
% See also: PSR_FT_JPSTH,  PSR_FT_CONVERT2FIELDTRIP

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

nUnits   = length(unitIDs);
nCombis  = 0.5 * (nUnits^2 - nUnits);
combis   = cell(nCombis,2);
itr      = 1;
for iUnit = 1:nUnits
    for jUnit = iUnit+1:nUnits
        combis(itr,:) = [{spikesFT.label{unitIDs(iUnit)}},{spikesFT.label{unitIDs(jUnit)}}];
        itr = itr + 1;
    end 
end

end