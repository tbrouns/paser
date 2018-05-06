function [spikes,freq] = psr_stimulus_window(spikes,freq,metadata,parameters)

Fs   = spikes.Fs;
Nmax = floor(Fs * metadata.duration) + 1;

% Extract stimulus windows

stimulusTimes = metadata.stimtimes{1};
stimulusType  = metadata.stimtimes{2};

if (~isempty(stimulusTimes))
    stimulusOnset = stimulusTimes(:,1) + parameters.analysis.stimoffset;
    stimulusOnset = stimulusOnset(~isnan(stimulusOnset));
end

nOnsets = length(stimulusOnset); % number of stimulus onsets
    
if (~isempty(stimulusOnset))
    
    %% Local field potential
    
    if (~isempty(freq))
        freq = psr_ft_convert2trials(freq,parameters,stimulusOnset,stimulusType); 
    end
    
    %% Spikes
        
    switch stimulusType
        case 'onset'
            tPre  = parameters.analysis.spk.t_win(1);
            tPost = parameters.analysis.spk.t_win(2);
            tPost = repmat(tPost,nOnsets,1);
        case 'interval'
            tPre  = 0;
            tPost = stimulusTimes(:,2) - stimulusOnset;
    end

    spikePoints = round(Fs * spikes.spiketimes) + 1;
    spikes.trials = false(nOnsets,length(spikes.spiketimes));
    spikes.info.trialonset = stimulusOnset + tPre;
    spikes.info.trialtime  = [repmat(tPre,length(tPost),1),tPost];
    
    if (~isempty(spikePoints))
        
        spikeSignal = zeros(Nmax,1);
        [spikePoints,uniqueIDs] = unique(spikePoints);
        spikeSignal(spikePoints) = 1;
                
        sPre  = round(Fs * tPre);  % pre- stimulus window
        sPost = round(Fs * tPost); % post-stimulus window
        
        stimulusPoints = round(Fs * stimulusOnset) + 1;   

        for iTrial = 1:nOnsets

            % In-line version of "psr_get_spike_ids" for speed-up
            
            sOnset = stimulusPoints(iTrial);
            I = sOnset + sPre : sOnset + sPost(iTrial);
            I = I(I >= 1 & I <= Nmax);
            I = I(spikeSignal(I) == 1);
            spikeSignal(I) = spikeSignal(I) + 1;
            spikeIDs = spikeSignal(spikeSignal > 0) == 2;
            spikeSignal(I) = spikeSignal(I) - 1; % Reset for next iteration
            spikeIDs = uniqueIDs(spikeIDs);

            spikes.trials(iTrial,spikeIDs) = true;
        end

        del = ~any(spikes.trials,1); % remove spikes outside of trials
        spikes = psr_sst_remove_spikes(spikes,del,'delete');
    
    end
end

end