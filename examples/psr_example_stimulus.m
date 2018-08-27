function stimTimes = psr_example_stimulus(cfg)

    loadPath   = cfg.loadpath;
    parameters = cfg.parameters;
    
    stimTimes(:,1) = psr_ms_detect_onset(loadPath,parameters);
    stimTimes(:,2) = cellstr('onset');

end