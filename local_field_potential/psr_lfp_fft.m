function fftfreq = psr_lfp_fft(freq,parameters)

% PSR_LFP_FFT - Fourier transform using FieldTrip's FT_FREQANALYSIS
%
% Syntax:  fftfreq = psr_lfp_fft(freq,parameters)
%
% Inputs:
%    freq       - FieldTrip LFP data structure (see README)
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    fftfreq - Output from FT_FREQANALYSIS
%
% See also: FT_FREQANALYSIS

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

freq = psr_ft_nan_removal(freq);

cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmfft';
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

if (~isempty_field(parameters,'parameters.analysis.fft.foi'));       cfg.foi       = parameters.analysis.fft.foi;       end
if (~isempty_field(parameters,'parameters.analysis.fft.taper'));     cfg.taper     = parameters.analysis.fft.taper;     end
if (~isempty_field(parameters,'parameters.analysis.fft.pad'));       cfg.pad       = parameters.analysis.fft.pad;       end
if (~isempty_field(parameters,'parameters.analysis.fft.tapsmofrq')); cfg.tapsmofrq = parameters.analysis.fft.tapsmofrq; end
if (~isempty_field(parameters,'parameters.analysis.fft.keepchans')); cfg.keepchans = parameters.analysis.fft.keepchans; end

% Temporariry remove some fields
[freq,~] = psr_remove_field(freq,'artifacts');
[freq,~] = psr_remove_field(freq,'missing');

% Fast Fourier Transform
try    fftfreq = ft_freqanalysis(cfg,freq);
catch; fftfreq = []; return;
end

% Output results
if (~cfg.keepchans) % Average over probe channels
    fftfreq.powspctrm = nanmean(fftfreq.powspctrm,2);
    fftfreq.label     = fftfreq.label{1};
end

fftfreq.freq      = single(fftfreq.freq);
fftfreq.powspctrm = single(fftfreq.powspctrm);
fftfreq           = orderfields(fftfreq);

end