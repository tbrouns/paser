function data = psr_bp_filter(data,cfg)
    [B,A] = butter(cfg.order,[cfg.lower cfg.upper]/(cfg.Fs/2),'bandpass');
    data  = filtfilt(B,A,data); % Zero-phase digital filtering
end