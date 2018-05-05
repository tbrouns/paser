function [bandpwr,missing] = psr_lfp_bpw(data,parameters)

Fr      = parameters.Fr;
fRange  = parameters.analysis.bpw.frange;
tRange  = [];
fNorm   = [];
bandpwr = [];
missing = [];
    
if (~psr_isempty_field(parameters,'parameters.analysis.bpw.trange')); tRange = parameters.analysis.bpw.trange; end
if (~psr_isempty_field(parameters,'parameters.analysis.bpw.fnorm'));  fNorm  = parameters.analysis.bpw.fnorm;  end

if (psr_isempty_field(data,'data.trial') || psr_isempty_field(data,'data.time')); return; end

% Set parameters

nFreqs  =  size(fRange,1);
nTrials =  size(data.trial,2);
missing = zeros(nTrials,1);
bandpwr = zeros(nTrials,nFreqs);

for iTrial = 1:nTrials
    
    T = data.time {iTrial};
    X = data.trial{iTrial}';
    
    if (~isempty(tRange))    
        i1 = find(T >= tRange(1),1);
        i2 = find(T >  tRange(2),1) - 1;
        if (isempty(i1)); i1 = 1;         end
        if (isempty(i2)); i2 = length(T); end
        X = X(i1:i2,:);
        T = T(i1:i2);
    end
        
    missing(iTrial) = sum(isnan(sum(X,2))) / size(X,1);
    
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
        
        for iFreq = 1:nFreqs
            bandpwr(iTrial,iFreq) = mean(bandpower(X,Fr,fRange(iFreq,:))) / norm; % Mean across channels
        end
        
    end
end

end