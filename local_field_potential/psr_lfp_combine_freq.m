function timefreq = psr_lfp_combine_freq(timefreq_array)

% PSR_LFP_COMBINE_FREQ - Combine multiple power spectrum structures into one
%
% Syntax:  timefreq = psr_lfp_combine_freq(timefreq_array)
%
% Inputs:
%    timefreq_array - Cell array of output structures from PSR_LFP_TFA
%
% Outputs:
%    timefreq - Single data structure containing the data from each element
%               of the input cell array
%
% See also: PSR_LFP_TFA

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

nBlocks = length(timefreq_array);
timefreq = timefreq_array{1};
if (isempty(timefreq)); return; end
nDims   = ndims(timefreq.powspctrm);

% Initialize
sz    = zeros(1,nDims);
sz(1) = length(timefreq_array);
for iBlock = 1:nBlocks
    for iDim = 2:nDims
        if (isfield(timefreq_array{iBlock},'powspctrm'))
            n = size(timefreq_array{iBlock}.powspctrm,iDim);
            if (n > sz(iDim)); sz(iDim) = n; end            
        end
    end
end
timefreq.powspctrm = NaN(sz);

% Insert data
itr = 1;
for iBlock = 1:nBlocks
    if (isfield(timefreq_array{iBlock},'powspctrm'))
        x = timefreq_array{iBlock}.powspctrm;
        sz1 = size(x,    1);
        sz2 = size(x,nDims);
        I = itr : itr + sz1 - 1;
        if (nDims == 4); timefreq.powspctrm(I,:,:,1:sz2) = x;
        else,            timefreq.powspctrm(I,:,  1:sz2) = x;
        end
        itr = itr + sz1;
    end
end

% Change time vector
if (isfield(timefreq,'time'))
    T  = timefreq.time;
    dt = mean(diff(mean(T,1)));
    n  = size(timefreq.powspctrm,nDims);
    t  = T(1):dt:(n-1)*dt+T(1);
    timefreq.time = t;
end

end