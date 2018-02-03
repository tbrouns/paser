function output = psr_lfp_bandpower(data,cfg)

Fr     = cfg.Fr;
fRange = cfg.frange;

if (isfield(cfg,'fnorm')); fNorm = cfg.fnorm;
else,                      fNorm = [];
end

% Set NaNs to zero

nTrials = size(data.trial,2);
for iTrial = 1:nTrials
    dataTrial = sum(data.trial{iTrial});
    data.trial{iTrial}(:,isnan(dataTrial)) = 0;
end

% Set parameters

nFreqs = size(fRange,1);
output = zeros(nFreqs,1);

X = data.trial{1}';

% Normalization

if (~isempty(fNorm)); norm = mean(bandpower(X,Fr,fNorm));
else,                 norm = 1;
end


for iFreq = 1:nFreqs
    output(iFreq) = mean(bandpower(X,Fr,fRange(iFreq,:))) / norm; % Mean across channels
end

end