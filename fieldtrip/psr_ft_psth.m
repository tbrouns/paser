function psth = psr_ft_psth(spikesFT,parameters,trialIDs)

cfg            = [];
cfg.keeptrials = 'yes';
if (nargin == 3); cfg.trials = trialIDs; end

if (~psr_isempty_field(parameters,'parameters.analysis.psth.binsize'));    cfg.binsize    = parameters.analysis.psth.binsize;    end
if (~psr_isempty_field(parameters,'parameters.analysis.psth.outputunit')); cfg.outputunit = parameters.analysis.psth.outputunit; end

psth = ft_spike_psth(cfg,spikesFT);
    
if (parameters.analysis.psth.crop); psth = cropPSTH(psth); end

end

function psth = cropPSTH(psth)

trials = psth.trial;
trials = sum(trials,1);
trials = sum(trials,2);
trials = trials(:);
trials = isnan(trials);

n = find(trials,1);

if (~isempty(n))
    
    psth.avg   = psth.avg  (:,  1:n);
    psth.var   = psth.var  (:,  1:n);
    psth.dof   = psth.dof  (:,  1:n);
    psth.time  = psth.time (:,  1:n);
    psth.trial = psth.trial(:,:,1:n);
    
    psth.cfg.latency = [psth.time(1) psth.time(end)];
    
end

end