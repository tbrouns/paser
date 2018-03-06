function stimTimes = psr_avp_stimulus(loadPath,parameters)

    k = psr_session_strcmp(parameters,'passive');
    loadPathActive  = loadPath(~k);
    loadPathPassive = loadPath( k);

    %% Camera onset times

    % For active trials
    stimTimes(~k,1) = psr_cam_detection(loadPathActive);
    stimTimes(~k,2) = cellstr('interval');

    %% Magnetic field artifacts [MFA]

    % For passive trials
    if (any(k))
        parameters.general.stims = parameters.general.stimuli{k};
        stimTimes(k,1) = psr_ms_detect_onset(loadPathPassive,parameters);
        stimTimes(k,2) = cellstr('onset');
    end

end