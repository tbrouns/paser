function [out_spikes,out_freq] = psr_stimulus_window(spikes,freq,stimtimes,parameters)

% PSR_STIMULUS_WINDOW - Extract data window around stimulus onsets.
%
% Syntax:  [out_spikes,out_freq] = psr_stimulus_window(spikes,freq,stimtimes,parameters)
%
% Inputs:
%    spikes     - See README
%    freq       - See README
%    stimtimes  - Two-element cell array, structured as followed:
%                 stimtimes{1}: stimulus timings [sec]
%                 stimtimes{2}: 'onset' or 'interval'
% 
%                 ## If stimtimes{2} = 'onset', 
%                    then stimtimes{1} should be a column vector containing
%                    each stimulus onset
% 
%                 ## If stimtimes{2} = 'interval',
%                    then stimtimes{1} should be a two-column matrix
%                    containing the stimulus onsets in the first column and
%                    the stimulus offsets in the second column
% 
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    out_spikes - Contains the "trials" field, which indicates the stimulus
%                 windows
%    out_freq   - LFP data aligned to stimulus onset

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Output
out_spikes = [];
out_freq   = [];

% Extract stimulus windows

stimulusTimes  = stimtimes{1};
stimulusType   = stimtimes{2};
stimulusOnset  = [];
stimulusOffset = [];
stimTimesTemp  = [];

if (~isempty(stimulusTimes))
    stimulusOnset  = stimulusTimes(:,1) + parameters.analysis.stimoffset;
    keep = ~isnan(stimulusOnset);
    stimulusOnset = stimulusOnset(keep);
    stimTimesTemp = stimulusOnset;
end 

if (~isempty(stimulusTimes) && size(stimulusTimes,2) >= 2)
    stimulusOffset = stimulusTimes(:,2); 
    stimulusOffset = stimulusOffset(keep);
    stimTimesTemp  = [stimulusOnset,stimulusOffset]; 
end

stimulusTimes = stimTimesTemp;
nOnsets = size(stimulusTimes,1); % number of stimulus onsets

if (~isempty(stimulusTimes))
    
    %% Local field potential
    if (~isempty(freq)); out_freq = psr_ft_convert2trials(freq,parameters,stimulusTimes,stimulusType); end
    
    %% Spikes
    if (~isempty(spikes))
        
        spikeTimes = spikes.spiketimes;
        
        % Initialize
        spikes.info.trialtime  = [];
        spikes.info.trialonset = [];
        spikes.info.stimdur    = [];
        
        switch stimulusType
            case 'onset'
                tPre  = parameters.analysis.spk.t_win(1);
                tPost = parameters.analysis.spk.t_win(2);
                tPost = repmat(tPost,nOnsets,1);
            case 'interval'
                tPre  = 0;
                tPost = stimulusOffset - stimulusOnset;
        end

        spikes.info.trialonset = stimulusOnset;
        spikes.info.trialtime  = [repmat(tPre,length(tPost),1),tPost];

        trials = false(nOnsets,length(spikeTimes));
        
        if (~isempty(spikeTimes))
            for iTrial = 1:nOnsets
                tOnset = stimulusOnset(iTrial);
                spikeIDs = spikeTimes >= tOnset + tPre & spikeTimes <= tOnset + tPost(iTrial);
                trials(iTrial,spikeIDs) = true;
            end
            keep = any(trials,1); % only keep spikes inside of trials
            spikes = psr_sst_remove_spikes(spikes,keep,'keep');
            spikes.trials = trials(:,keep);
        end
        
        out_spikes = spikes;
    end
end

end