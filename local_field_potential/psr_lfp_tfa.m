function timefreq = psr_lfp_tfa(freq,parameters)

% PSR_LFP_TFA - Time-frequency analysis using FieldTrip's FT_FREQANALYSIS
%
% Syntax:  timefreq = psr_lfp_tfa(freq,parameters)
%
% Inputs:
%    freq       - FieldTrip LFP data structure (see README)
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    timefreq - Output from FT_FREQANALYSIS
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
[freq,~]       = psr_remove_field(freq,'artifacts');
[freq,missing] = psr_remove_field(freq,'missing');

% Time-frequency analysis
try timefreq = ft_freqanalysis(cfg,freq); % FieldTrip function
catch ME
    timefreq = [];
    str = ME.message;
    psr_show_warning({str});
    return;
end

% Set artifacts to NaN
nTime = length(timefreq.time);
nChan =   size(timefreq.powspctrm,2);
for iTrial = 1:length(missing)
    for iChan = 1:nChan
        I  = missing{iTrial}(iChan,:);
        N  = size(I,2);
        I  = find(I);
        i  = unique(ceil((nTime/N)*I));        
        timefreq.powspctrm(iTrial,iChan,:,i) = NaN;
    end
end
   
% Output results
if (~cfg.keepchans) % Average over probe channels
    timefreq.powspctrm = nanmean(timefreq.powspctrm,2);
    timefreq.label     = timefreq.label{1};
end

timefreq.powspctrm = single(timefreq.powspctrm); % To avoid memory issues
timefreq = orderfields(timefreq);

end