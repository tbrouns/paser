function bandpwr = psr_lfp_bpw(freq,parameters)

% PSR_LFP_BPW - Bandpower calculation from time-series data
%
% Syntax:  bandpwr = psr_lfp_bpw(freq,parameters)
%
% Inputs:
%    freq       - FieldTrip LFP data structure
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    bandpwr - Array containing bandpowers, with shape:
%              [Number of trials x Number of frequency bands]
%
% Dependencies: Signal Processing Toolbox
%
% See also: BANDPOWER

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

Fr      = parameters.Fr;
fRange  = parameters.analysis.bpw.frange;
tRange  = [];
fNorm   = [];
bandpwr = [];
    
if (~isempty_field(parameters,'parameters.analysis.bpw.trange')); tRange = parameters.analysis.bpw.trange; end
if (~isempty_field(parameters,'parameters.analysis.bpw.fnorm'));  fNorm  = parameters.analysis.bpw.fnorm;  end

if (isempty_field(freq,'freq.trial') || isempty_field(freq,'freq.time')); return; end

% Set parameters

nBands  =  size(fRange,1);
nTrials =  size(freq.trial,2);
bandpwr = zeros(nTrials,nBands);

for iTrial = 1:nTrials
    
    T = freq.time {iTrial};
    X = freq.trial{iTrial}';
    
    if (~isempty(tRange))    
        i1 = find(T >= tRange(1),1);
        i2 = find(T >  tRange(2),1) - 1;
        if (isempty(i1)); i1 = 1;         end
        if (isempty(i2)); i2 = length(T); end
        X = X(i1:i2,:);
        T = T(i1:i2);
    end
            
    % Remove NaNs
    
    temp = [];
    temp.trial = {X'};
    temp.time  = {T};
    temp = psr_ft_nan_removal(temp);
    X = temp.trial{1}';
    
    if (~isempty(X))
        
        if (~isempty(fNorm)); norm = mean(bandpower(X,Fr,fNorm)); % Normalization
        else,                 norm = 1;
        end
        
        for iFreq = 1:nBands
            bandpwr(iTrial,iFreq) = mean(bandpower(X,Fr,fRange(iFreq,:))) / norm; % Mean across channels
        end
        
    end
end

end