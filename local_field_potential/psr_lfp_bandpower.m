function [bandpwr,missing] = psr_lfp_bandpower(data,parameters)

Fr     =  parameters.Fr;
fRange =  parameters.analysis.bpw.frange;

str = 'parameters.analysis.bpw.trange';
if (~psr_isempty_field(parameters,str)); tRange = parameters.analysis.bpw.trange;
else,                                    tRange = [];
end

str = 'parameters.analysis.bpw.fnorm';
if (~psr_isempty_field(parameters,str)); fNorm = parameters.analysis.bpw.fnorm;
else,                                    fNorm = [];
end

% Set parameters

nFreqs  =  size(fRange,1);
nTrials =  size(data.trial,2);
missing = zeros(nTrials,1);
bandpwr = zeros(nTrials,nFreqs);

for iTrial = 1:nTrials
    
    T = data.time{iTrial};
    X = data.trial{iTrial}';
    
    if (~isempty(tRange))    
        i1 = find(T >= tRange(1),1);
        i2 = find(T >  tRange(2),1) - 1;
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
    
    % Normalization

    if (~isempty(fNorm)); norm = mean(bandpower(X,Fr,fNorm));
    else,                 norm = 1;
    end

    for iFreq = 1:nFreqs
        bandpwr(iTrial,iFreq) = mean(bandpower(X,Fr,fRange(iFreq,:))) / norm; % Mean across channels
    end
end

end