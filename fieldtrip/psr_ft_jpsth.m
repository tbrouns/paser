function [jpsth,jpsthShuff] = psr_ft_jpsth(psth,parameters,clustIDs)

% Wrapper function
    
cfg        = [];
cfg.method = 'jpsth';

if (nargin == 3); cfg.channelcmb = psr_ft_combis(psth,clustIDs); end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.keeptrials')); cfg.keeptrials    = parameters.analysis.jpsth.keeptrials; end
if (~psr_isempty_field(parameters,'parameters.analysis.jpsth.normalize'));  cfg.normalization = parameters.analysis.jpsth.normalize;  end

jpsth = psr_ft_spike_jpsth(cfg,psth);

% Subtract the shift predictor
jpsthShuff = [];
if (parameters.analysis.jpsth.shuffle)
    cfg.method = 'shiftpredictor';
    jpsthShuff  = psr_ft_spike_jpsth(cfg,psth);
    jpsth.jpsth = jpsth.jpsth - jpsthShuff.shiftpredictor;
end

end