function output = psr_ft_cat_trials(spikesFT,twin)

% spikesFT: fieldTrip structure containing spike info
% twin: window of each trial that you wish to concatenate

if (nargin < 2); twin = []; end

nUnits  = length(spikesFT.label);
nTrials = size(spikesFT.trialtime,1);

spiketimes = cell(1,nUnits);
trialIDs   = cell(1,nUnits);
trialonset = [0;cumsum(diff(spikesFT.trialtime,1,2))];

for iUnit = 1:nUnits
    spiketimesUnit = cell(1,nTrials);
    trialIDsUnit   = cell(1,nTrials);
    for iTrial = 1:nTrials
        spikeIDs       = spikesFT.trial{iUnit} == iTrial;
        spiketimesTemp = spikesFT.time{iUnit}(spikeIDs);
        if (~isempty(twin)); spiketimesTemp = spiketimesTemp(spiketimesTemp >= twin(1) & spiketimesTemp <= twin(2)); end
        spiketimesTemp = spiketimesTemp + trialonset(iTrial) - spikesFT.trialtime(iTrial,1);
        spiketimesUnit{iTrial} = spiketimesTemp;
        trialIDsUnit  {iTrial} = iTrial * ones(size(spiketimesTemp));
    end
    spiketimes{iUnit} = cell2mat(spiketimesUnit);
    trialIDs  {iUnit} = cell2mat(trialIDsUnit);
end

tobs = trialonset(end);

output       = [];
output.time  = spiketimes;
output.trial = trialIDs;
output.tobs  = tobs;
output.units = spikesFT.label;

end