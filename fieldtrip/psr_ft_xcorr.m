function XCorr = psr_ft_xcorr(spikesFT,parameters)

cfg            = [];
cfg.method     = 'xcorr';
cfg.keeptrials = 'yes';

if (~isempty_field(parameters,'parameters.analysis.xcorr.binsize'));    cfg.binsize    = parameters.analysis.xcorr.binsize;    end
if (~isempty_field(parameters,'parameters.analysis.xcorr.debias'));     cfg.debias     = parameters.analysis.xcorr.debias;     end
if (~isempty_field(parameters,'parameters.analysis.xcorr.maxlag'));     cfg.maxlag     = parameters.analysis.xcorr.maxlag;     end
if (~isempty_field(parameters,'parameters.analysis.xcorr.outputunit')); cfg.outputunit = parameters.analysis.xcorr.outputunit; end
    
XCorr = psr_ft_spike_xcorr(cfg,spikesFT);

% Subtract the shift predictor
if (parameters.analysis.xcorr.shuffle) 
    cfg.method  = 'shiftpredictor';
    Xshuff      = psr_ft_spike_xcorr(cfg,spikesFT);
    XCorr.xcorr = XCorr.xcorr - Xshuff.shiftpredictor;
end

end