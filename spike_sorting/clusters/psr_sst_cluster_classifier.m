function spikes = psr_sst_cluster_classifier(spikes,parameters)

%% Neural network classification

rpvs   = {spikes.clusters.metrics.rpv};
pAuc   = {spikes.clusters.metrics.p_auc};
ampRel = {spikes.clusters.metrics.amp_rel};
corrG  = {spikes.clusters.metrics.corrG};
corrP  = {spikes.clusters.metrics.corrP};
subs   = {spikes.clusters.metrics.sub};
xclag  = {spikes.clusters.metrics.xc_lag};
snrs   = {spikes.clusters.metrics.snr};

I = ...
    ~cellfun(@isempty,rpvs)   & ...
    ~cellfun(@isempty,pAuc)   & ...
    ~cellfun(@isempty,ampRel) & ...
    ~cellfun(@isempty,corrG)  & ...
    ~cellfun(@isempty,corrP)  & ...
    ~cellfun(@isempty,subs)   & ...
    ~cellfun(@isempty,xclag)  & ...
    ~cellfun(@isempty,snrs);

I = find(I);

rpvs   = cell2mat(  rpvs(I));
pAuc   = cell2mat(  pAuc(I));
ampRel = cell2mat(ampRel(I));
corrG  = cell2mat( corrG(I));
corrP  = cell2mat( corrP(I));
subs   = cell2mat(  subs(I)); 
xclag  = cell2mat( xclag(I));
snrs   = cell2mat(  snrs(I));

X = [rpvs;pAuc;ampRel;corrG;corrP;subs;xclag;snrs];
load('clusterClassifier.mat'); % load neural network
labels = net(X);

flag = ones(size(I));

%% Initialize types to zero
for iClust = 1:size(spikes.clusters.metrics,2); spikes.clusters.metrics(iClust).type = 0; end

%% Hard thresholds
M = {spikes.clusters.metrics.amp_rel};  flag(cell2mat(M(I)) < 1.0)                              = 0; % mean sub-threshold amplitude
M = {spikes.clusters.metrics.amp};      flag(cell2mat(M(I)) > parameters.cluster.max_amplitude) = 0; % absolute amplitude
M = {spikes.clusters.metrics.p2p};      flag(cell2mat(M(I)) > parameters.cluster.max_p2p)       = 0; % absolute peak-to-peak
M = {spikes.clusters.metrics.sub};      flag(cell2mat(M(I)) > parameters.cluster.max_sub)       = 0; % fraction of sub-threshold spikes
M = {spikes.clusters.metrics.rpv};      flag(cell2mat(M(I)) > parameters.cluster.max_rpv)       = 0; % fraction of refractory period violations
M = {spikes.clusters.metrics.nspikes};  flag(cell2mat(M(I)) < parameters.cluster.min_spikes)    = 0; % small spike number
M = {spikes.clusters.metrics.frate};    flag(cell2mat(M(I)) < parameters.cluster.min_frate)     = 0; % firing rate
M = {spikes.clusters.metrics.Lratio};   flag(cell2mat(M(I)) > parameters.cluster.max_lratio)    = 0; % l-ratio
M = {spikes.clusters.metrics.IsoDis};   flag(cell2mat(M(I)) < parameters.cluster.min_isodist)   = 0; % isolation distance

%% Set cluster classes
nClusts = length(flag);
labels  = flag .* labels;
for iflag = 1:nClusts
    spikes.clusters.metrics(I(iflag)).type = labels(iflag);
end

end