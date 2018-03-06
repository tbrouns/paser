function output = psr_lfp_tfa(data,parameters)

% Data preparation for FT_FREQANALYSIS

data = psr_ft_nan_removal(data);

cfg            = [];
cfg.Fs         = parameters.Fr;
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

if (~psr_isempty_field(parameters,'parameters.analysis.tfa.taper')); cfg.taper = parameters.analysis.tfa.taper; end
if (~psr_isempty_field(parameters,'parameters.analysis.tfa.pad'));   cfg.pad   = parameters.analysis.tfa.pad;   end

fLower  = parameters.analysis.tfa.freq.lower;
fUpper  = parameters.analysis.tfa.freq.upper;
fStep   = parameters.analysis.tfa.freq.step;
cfg.foi = fLower:fStep:fUpper;
nfoi    = length(cfg.foi);

if (isfield(parameters.analysis.tfa,'ncycles') && ~isempty(parameters.analysis.tfa.ncycles))
    cfg.t_ftimwin = parameters.analysis.tfa.ncycles ./ cfg.foi;
else
    cfg.t_ftimwin = parameters.analysis.tfa.twin * ones(1,nfoi);
end

if (size(data.time,2) > 1)
    onset  = parameters.lfp.trial.onset;
    offset = parameters.lfp.trial.offset;
else % Single continuous trial
    onset  = data.time{1}(1);
    offset = data.time{1}(end);
end

cfg.toi = onset:parameters.analysis.tfa.time.step:offset;

switch parameters.analysis.tfa.smtype % Smoothing
    case 'freq';  cfg.tapsmofrq = parameters.analysis.tfa.smfreq  * ones(1,nfoi);
    case 'ratio'; cfg.tapsmofrq = parameters.analysis.tfa.smratio * cfg.foi;
end

% Temporariry remove some fields
[data,~]       = psr_remove_field(data,'artifacts');
[data,missing] = psr_remove_field(data,'missing');

% Time-frequency analysis
output = ft_freqanalysis(cfg,data);
    
% Set artifacts to NaN
n = length(output.time);
for iTrial = 1:length(missing)
    T  = data.time{iTrial};
    I  =   missing{iTrial};
    i1 = find(T >= onset,1);
    i2 = find(T > offset,1) - 1;
    I  =  I(i1:i2);
    N  = length(I);
    I  =   find(I);
    i  = unique(ceil((n/N)*I));        
    output.powspctrm(iTrial,:,:,i) = NaN;
end

% Output results
if (~parameters.analysis.tfa.keepchans) % Average over probe channels
    output.powspctrm = nanmean(output.powspctrm,2);
    output.label     = output.label{1};
end

output.freq      = single(output.freq);
output.powspctrm = single(output.powspctrm);
output.time      = single(output.time); 
output           = orderfields(output);

end