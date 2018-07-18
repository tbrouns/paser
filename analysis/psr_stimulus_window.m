function [out_spikes,out_freq] = psr_stimulus_window(spikes,freq,metadata,parameters)

% Output
out_spikes = [];
out_freq   = [];

% Extract stimulus windows

stimulusTimes  = metadata.stimtimes{1};
stimulusType   = metadata.stimtimes{2};
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