function freq = psr_lfp_timefreq(data,parameters)

onset  = parameters.lfp.trial_onset  - parameters.lfp.trial_padding;
offset = parameters.lfp.trial_offset + parameters.lfp.trial_padding;

cfg            = [];
cfg.method     = parameters.lfp.method;
cfg.taper      = parameters.lfp.taper;
cfg.output     = 'pow';
cfg.pad        = parameters.lfp.pad;
if (strcmp(cfg.method,'mtmfft'))
    cfg.foilim     = [parameters.lfp.freq_lower parameters.lfp.freq_upper];
else
    cfg.foi        = parameters.lfp.freq_lower:parameters.lfp.freq_step:parameters.lfp.freq_upper;
    cfg.toi        = onset:parameters.lfp.time_step:offset;
    cfg.t_ftimwin  = parameters.lfp.ncycles ./ cfg.foi;  % 5 cycles per (sliding) time window
    cfg.t_ftimwin(cfg.t_ftimwin > parameters.lfp.trial_padding) = parameters.lfp.trial_padding;
end
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
freq = ft_freqanalysis(cfg, data);

% Remove trials with excessive data gaps

psr_parameter_config; % TEMP

nTrials = size(freq.powspctrm,1);
del = false(nTrials,1);
for iTrial = 1:nTrials
    pow = squeeze(freq.powspctrm(iTrial,:,:,:));
    Ntot = sum(sum(sum(~isnan(pow))));
    N    = sum(sum(sum(pow < 0.001)));
    if (N / Ntot > parameters.lfp.miss_thresh)
        del(iTrial) = true;
    end
end

freq.powspctrm(del,:,:,:) = []; % Delete
freq.pow       = freq.powspctrm; % TEMP
freq.trialIDs  = del;
freq.std       = squeeze( nanstd(freq.powspctrm,[],1)); % Standard deviation for every bin
freq.powspctrm = squeeze(nanmean(freq.powspctrm,1)); % Average over trials

end