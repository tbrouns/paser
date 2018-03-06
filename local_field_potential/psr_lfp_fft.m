function output = psr_lfp_fft(data,parameters)

% Data preparation for FT_FREQANALYSIS

data = psr_ft_nan_removal(data);

cfg            = [];
cfg.Fs         = parameters.Fr;
cfg.output     = 'pow';
cfg.method     = 'mtmfft';
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

if (~psr_isempty_field(parameters,'parameters.analysis.fft.taper'));     cfg.taper     = parameters.analysis.fft.taper;     end
if (~psr_isempty_field(parameters,'parameters.analysis.fft.pad'));       cfg.pad       = parameters.analysis.fft.pad;       end
if (~psr_isempty_field(parameters,'parameters.analysis.fft.tapsmofrq')); cfg.tapsmofrq = parameters.analysis.fft.tapsmofrq; end

fLower  = parameters.analysis.fft.freq.lower;
fUpper  = parameters.analysis.fft.freq.upper;
fStep   = parameters.analysis.fft.freq.step;
cfg.foi = fLower:fStep:fUpper;

% Temporariry remove some fields
[data,~] = psr_remove_field(data,'artifacts');
[data,~] = psr_remove_field(data,'missing');

% Fast Fourier Transform
output = ft_freqanalysis(cfg,data);
    
% Output results
if (~parameters.analysis.fft.keepchans) % Average over probe channels
    output.powspctrm = nanmean(output.powspctrm,2);
    output.label     = output.label{1};
end

output.freq      = single(output.freq);
output.powspctrm = single(output.powspctrm);
output           = orderfields(output);

end