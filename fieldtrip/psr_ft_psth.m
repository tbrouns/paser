function psth = psr_ft_psth(spikesFT,parameters,trialIDs)

% PSR_FT_PSTH - Wrapper function for FieldTrip's FT_SPIKE_PSTH
% This function calls the FieldTrip function: FT_SPIKE_PSTH, which
% calculates peri-stimulus time histogram (PSTH) metrics
% 
% Syntax:  psth = psr_ft_psth(spikesFT,parameters,trialIDs)
%
% Inputs:
%    spikesFT   - FieldTrip spike structure
%    parameters - See PSR_PARAMETERS_ANALYSIS
%    trialIDs   - Which trials to include in the PSTH calculation
%
% Outputs:
%    psth - Output from FT_SPIKE_PSTH
%
% See also: FT_SPIKE_PSTH, PSR_FT_CONVERT2FIELDTRIP

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Initialize
psth           = [];
cfg            = [];
cfg.keeptrials = 'yes';
if (nargin == 3); cfg.trials = trialIDs; end
if (isempty(spikesFT)); return; end
if (~isempty_field(parameters,'parameters.analysis.psth.binsize'));    cfg.binsize    = parameters.analysis.psth.binsize;    end
if (~isempty_field(parameters,'parameters.analysis.psth.outputunit')); cfg.outputunit = parameters.analysis.psth.outputunit; end
if (isfield(cfg,'binsize') && any(spikesFT.trialtime(:,2) < cfg.binsize)); return; end

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