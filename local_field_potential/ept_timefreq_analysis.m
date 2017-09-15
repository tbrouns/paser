function freq = ept_timefreq_analysis(data,parameters)
    
    tlength = parameters.lfp.trial_length - parameters.lfp.base_onset; 

    cfg            = [];
    cfg.method     = parameters.lfp.method;
    cfg.taper      = parameters.lfp.taper;
    cfg.output     = 'pow';
    cfg.pad        = parameters.lfp.pad;
    if (strcmp(cfg.method,'mtmfft'))
        cfg.foilim     = [parameters.lfp.freq_lower parameters.lfp.freq_upper];
    else
        cfg.foi        = parameters.lfp.freq_lower:parameters.lfp.freq_step:parameters.lfp.freq_upper;
        cfg.toi        = parameters.lfp.base_onset:parameters.lfp.time_step:parameters.lfp.trial_length;
%         cfg.t_ftimwin  = ones(size(cfg.foi)).*0.5;
        cfg.t_ftimwin  = 5./cfg.foi;  % 5 cycles per (sliding) time window 
        cfg.t_ftimwin(cfg.t_ftimwin > tlength) = tlength;
    end
    cfg.channel    = 'all';
    cfg.trials     = 'all';
    cfg.keeptrials = 'no';
    cfg.keeptapers = 'no';
    freq           = ft_freqanalysis(cfg, data);

end