function output = psr_lfp_tfa(data,parameters)

% Data preparation for FT_FREQANALYSIS

data = psr_ft_nan_removal(data);

cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

if (~isempty_field(parameters,'parameters.analysis.tfa.tapsmofrq')); cfg.tapsmofrq = parameters.analysis.tfa.tapsmofrq; end
if (~isempty_field(parameters,'parameters.analysis.tfa.toi'));       cfg.toi       = parameters.analysis.tfa.toi;       end
if (~isempty_field(parameters,'parameters.analysis.tfa.foi'));       cfg.foi       = parameters.analysis.tfa.foi;       end
if (~isempty_field(parameters,'parameters.analysis.tfa.taper'));     cfg.taper     = parameters.analysis.tfa.taper;     end
if (~isempty_field(parameters,'parameters.analysis.tfa.pad'));       cfg.pad       = parameters.analysis.tfa.pad;       end
if (~isempty_field(parameters,'parameters.analysis.tfa.keepchans')); cfg.keepchans = parameters.analysis.tfa.keepchans; end
if (~isempty_field(parameters,'parameters.analysis.tfa.ncycles'));   cfg.t_ftimwin = parameters.analysis.tfa.ncycles ./ cfg.foi;
else,                                                                cfg.t_ftimwin = parameters.analysis.tfa.t_ftimwin * ones(size(cfg.foi));
end

% Temporariry remove some fields
[data,~]       = psr_remove_field(data,'artifacts');
[data,missing] = psr_remove_field(data,'missing');

% Time-frequency analysis
try output = ft_freqanalysis(cfg,data); % FieldTrip function
catch ME
    output = [];
    str = ME.message;
    psr_show_warning({str});
    return;
end

% Set artifacts to NaN
nTime = length(output.time);
nChan =   size(output.powspctrm,2);
for iTrial = 1:length(missing)
    for iChan = 1:nChan
        I  = missing{iTrial}(iChan,:);
        N  = size(I,2);
        I  = find(I);
        i  = unique(ceil((nTime/N)*I));        
        output.powspctrm(iTrial,iChan,:,i) = NaN;
    end
end
   
% Output results
if (~cfg.keepchans) % Average over probe channels
    output.powspctrm = nanmean(output.powspctrm,2);
    output.label     = output.label{1};
end

output.powspctrm = single(output.powspctrm); % To avoid memory issues
output = orderfields(output);

end