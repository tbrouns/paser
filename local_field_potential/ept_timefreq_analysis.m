function freq = ept_timefreq_analysis(data,parameters,stim)

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
    cfg.keeptrials = 'no';
    cfg.keeptapers = 'no';
    freq           = ft_freqanalysis(cfg, data);

    if (stim == 0) % average over whole time window
        pow = freq.powspctrm;
        pow = nanmean(pow,3);
        freq.powspctrm = repmat(pow,1,1,size(freq.time,2));
    end
end