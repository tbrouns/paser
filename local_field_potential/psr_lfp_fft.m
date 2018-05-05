function output = psr_lfp_fft(data,parameters)

% Data preparation for FT_FREQANALYSIS

data = psr_ft_nan_removal(data);

cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmfft';
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

if (~psr_isempty_field(parameters,'parameters.analysis.fft.foi'));       cfg.foi       = parameters.analysis.fft.foi;       end
if (~psr_isempty_field(parameters,'parameters.analysis.fft.taper'));     cfg.taper     = parameters.analysis.fft.taper;     end
if (~psr_isempty_field(parameters,'parameters.analysis.fft.pad'));       cfg.pad       = parameters.analysis.fft.pad;       end
if (~psr_isempty_field(parameters,'parameters.analysis.fft.tapsmofrq')); cfg.tapsmofrq = parameters.analysis.fft.tapsmofrq; end
if (~psr_isempty_field(parameters,'parameters.analysis.fft.keepchans')); cfg.keepchans = parameters.analysis.fft.keepchans; end

% Temporariry remove some fields
[data,~] = psr_remove_field(data,'artifacts');
[data,~] = psr_remove_field(data,'missing');

% Fast Fourier Transform
try    output = ft_freqanalysis(cfg,data);
catch; output = []; return;
end

% Output results
if (~cfg.keepchans) % Average over probe channels
    output.powspctrm = nanmean(output.powspctrm,2);
    output.label     = output.label{1};
end

output.freq      = single(output.freq);
output.powspctrm = single(output.powspctrm);
output           = orderfields(output);

end