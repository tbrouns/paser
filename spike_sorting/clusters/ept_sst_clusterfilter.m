function spikes = ept_sst_clusterfilter(spikes,parameters)

flag = true(size([spikes.clusters.vars.id]));
type = zeros(size(flag));

flag([spikes.clusters.vars.nspikes]  < parameters.cluster.min_spikes) = 0; % small spike number
flag([spikes.clusters.vars.amp_rel]  < 1.0)                           = 0; % mean sub-threshold amplitude
    
type(type == 0 & ~flag) = 1; % ignore these clusters

flag([spikes.clusters.vars.artifact] > parameters.cluster.max_artifact)  = 0; % overlap with LFP artifact regions
flag([spikes.clusters.vars.amp]      > parameters.cluster.max_amplitude) = 0; % absolute amplitude

type(type == 0 & ~flag) = 2; % artifacts

flag([spikes.clusters.vars.rpv]      > parameters.cluster.max_rpv)       = 0; % fraction of refractory period violations
flag([spikes.clusters.vars.missing]  > parameters.cluster.max_missing)   = 0; % fraction of missing spikes
flag([spikes.clusters.vars.xc_lag]   > parameters.cluster.max_lag)       = 0; % lag between possible intertwined signals

type(type == 0 & ~flag) = 3; % multi-unit cluster

flag([spikes.clusters.vars.Lratio]   > parameters.cluster.max_lratio)    = 0; % l-ratio
flag([spikes.clusters.vars.IsoDis]   < parameters.cluster.min_isodist)   = 0; % isolation distance

type(type == 0 & ~flag) = 4; % partial single-unit cluster
type(type == 0)         = 5; % single cluster

nflags = length(flag);

for iflag = 1:nflags 
    spikes.clusters.vars(iflag).type = type(iflag);
end

end