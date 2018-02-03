function output = psr_lfp_timefreq(data,parameters,method)

% For details on the following parameters, see FT_FREQANALYSIS help.

if (nargin < 3); method = []; end

data = psr_ft_nan_removal(data);

params = parameters.analysis.tfa;

% Time-frequency analysis

cfg           = [];
cfg.Fs        = parameters.Fr;
cfg.output    = 'pow';

if (isfield(params,'method')); cfg.method = params.method; end
if (isfield(params,'taper'));  cfg.taper  = params.taper;  end
if (isfield(params,'pad'));    cfg.pad    = params.pad;    end

switch cfg.method
    
    case 'mtmfft'
        
        cfg.foi = params.freq.lower:params.freq.step:params.freq.upper;
        cfg.tapsmofrq = params.smfreq; % Has to be scalar
        
    case 'mtmconvol'
        
        cfg.foi = params.freq.lower:params.freq.step:params.freq.upper;
        nfoi = length(cfg.foi);
        
        if (isfield(params,'ncycles') && ~isempty(params.ncycles))
            cfg.t_ftimwin = params.ncycles ./ cfg.foi;
        else
            cfg.t_ftimwin = params.twin * ones(1,nfoi);
        end
        
        if (strcmp(method,'stimulus'))
            cfg.t_ftimwin(cfg.t_ftimwin > parameters.lfp.trial.padding) = parameters.lfp.trial.padding;
            onset  = parameters.lfp.trial.onset  - parameters.lfp.trial.padding;
            offset = parameters.lfp.trial.offset + parameters.lfp.trial.padding;
        else % Single continuous trial
            onset  = data.time{1}(1);
            offset = data.time{1}(end);
        end
        
        cfg.toi = onset:params.time.step:offset;
        
        switch params.smtype
            case 'freq';  cfg.tapsmofrq = params.smfreq  * ones(1,nfoi);
            case 'ratio'; cfg.tapsmofrq = params.smratio * cfg.foi;
        end
end

cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';

% Temporariry remove artifacts field
artifacts = data.artifacts;
data = rmfield(data,'artifacts');

% Frequency analysis
output = ft_freqanalysis(cfg, data);

% Set artifacts to NaN
if isfield(output,'time')
    output = artifact_removal(output,artifacts,cfg.Fs);
end

% Output results
if (~params.keepchans) % Average over probe channels
    output.powspctrm = mean(output.powspctrm,2);
end

output.artifacts = artifacts;
output.freq      = single(output.freq);
output.powspctrm = single(output.powspctrm);
if (isfield(output,'time')); output.time = single(output.time); end
output = orderfields(output);

end

function data = artifact_removal(data,artifacts,Fs)

artifacts = (artifacts - 1) / Fs;

% Ignore artifacts that are shorter than time bin size
dt = artifacts(:,2) - artifacts(:,1);
dT = mean(diff(data.time));
del = dt < dT;
artifacts(del,:) = [];

% Set artifacts to NaN
T = data.time(end);
artifacts = round(length(data.time) * artifacts / T) + 1;
for i = 1:size(artifacts,1)
    t = artifacts(i,1):artifacts(i,2);
    data.powspctrm(:,:,:,t) = NaN;
end

end
