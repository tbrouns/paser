function output = psr_ft_crop_trials(spikesFT,twin)

% PSR_FT_CROP_TRIALS - Crop all trials in FieldTrip spike structure
% With this function you can reduce the size of every trial in a FieldTrip
% spike data structure 
% 
% Syntax:  output = psr_ft_crop_trials(spikesFT,twin)
%
% Inputs:
%    spikesFT - A FieldTrip spike data structure
% 
%    twin     - Two-element vector indicating start and end point [sec] of
%               the new trial window
%
% Outputs:
%    output - Structure containing fields:
%             "time"  : Cell array with each cell containing a vector of
%                       spike times for the corresponding unit
%             "trial" : Cell array with each cell containing a vector of
%                       trial labels for the spike timestamps given in the
%                       "time" field
%             "tobs"  : Total observational period
%             "units" : Label for each unit
% 
% See also: PSR_FT_CONVERT2FIELDTRIP

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (nargin < 2); twin = []; end

nUnits  = length(spikesFT.label);
nTrials = size(spikesFT.trialtime,1);

spiketimes = cell(1,nUnits);
trialIDs   = cell(1,nUnits);
trialonset = [0;cumsum(diff(spikesFT.trialtime,1,2))];

for iUnit = 1:nUnits
    unit_spiketimes = cell(1,nTrials);
    unit_trialIDs   = cell(1,nTrials);
    for iTrial = 1:nTrials
        spikeIDs       = spikesFT.trial{iUnit} == iTrial; % Grab the spikes for the current trial
        spiketimesTemp = spikesFT.time{iUnit}(spikeIDs);
        if (~isempty(twin)); spiketimesTemp = spiketimesTemp(spiketimesTemp >= twin(1) & spiketimesTemp <= twin(2)); end % Crop the trial
        spiketimesTemp = spiketimesTemp + trialonset(iTrial) - spikesFT.trialtime(iTrial,1); % Offset spike times for new trial window
        unit_spiketimes{iTrial} = spiketimesTemp;
        unit_trialIDs  {iTrial} = iTrial * ones(size(spiketimesTemp));
    end
    spiketimes{iUnit} = cell2mat(unit_spiketimes);
    trialIDs  {iUnit} = cell2mat(unit_trialIDs);
end

tobs = trialonset(end);

output       = [];
output.time  = spiketimes;     % Spike times for each unit
output.trial = trialIDs;       % Trial IDs   for each unit
output.tobs  = tobs;           % Total duration of the session
output.units = spikesFT.label; % All unit IDs

end