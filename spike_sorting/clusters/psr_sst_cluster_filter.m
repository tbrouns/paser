function spikes = psr_sst_cluster_filter(spikes,parameters)

% Integers to classify cluster so it easier to extract clusters above
% certain 'quality' level
%
% Types:
% 1: Noise
% 2: Artifact
% 3: Multi-unit
% 4: Partial single-unit
% 5: Full single-unit

flag = true(size([spikes.clusters.vars.id]));
type = zeros(size(flag));

flag([spikes.clusters.vars.nspikes]  < parameters.cluster.min_spikes) = 0; % small spike number
flag([spikes.clusters.vars.amp_rel]  < 1.0)                           = 0; % mean sub-threshold amplitude
    
type(type == 0 & ~flag) = 1; % Noise

flag([spikes.clusters.vars.artifact] > parameters.cluster.max_artifact)  = 0; % overlap with LFP artifact regions
flag([spikes.clusters.vars.amp]      > parameters.cluster.max_amplitude) = 0; % absolute amplitude
flag([spikes.clusters.vars.corrG]    > parameters.cluster.max_corr)      = 0; % correlation across probes
flag([spikes.clusters.vars.p_auc]    < parameters.cluster.min_pauc)      = 0; % area under curve for Poisson distribution
flag([spikes.clusters.vars.frate]    < parameters.cluster.min_frate)     = 0; % firing rate

type(type == 0 & ~flag) = 2; % Artifacts

flag([spikes.clusters.vars.rpv]      > parameters.cluster.max_rpv)       = 0; % fraction of refractory period violations
flag([spikes.clusters.vars.missing]  > parameters.cluster.max_missing)   = 0; % fraction of missing spikes
flag([spikes.clusters.vars.xc_lag]   > parameters.cluster.max_lag)       = 0; % lag between possible intertwined signals

type(type == 0 & ~flag) = 3; % Multi-unit cluster

flag([spikes.clusters.vars.Lratio]   > parameters.cluster.max_lratio)    = 0; % l-ratio
flag([spikes.clusters.vars.IsoDis]   < parameters.cluster.min_isodist)   = 0; % isolation distance

type(type == 0 & ~flag) = 4; % Partial single-unit cluster
type(type == 0)         = 5; % Single-unit cluster

nflags = length(flag);

for iflag = 1:nflags 
    spikes.clusters.vars(iflag).type = type(iflag);
end

end