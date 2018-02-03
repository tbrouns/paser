function spikes = psr_sst_cluster_classifier(spikes)

%% Neural network classification

rpvs   = {spikes.clusters.metrics.rpv};
pAuc   = {spikes.clusters.metrics.cAuc};
ampRel = {spikes.clusters.metrics.ampRel};
corrG  = {spikes.clusters.metrics.corrG};
corrP  = {spikes.clusters.metrics.corrP};
subs   = {spikes.clusters.metrics.sub};
xcLag  = {spikes.clusters.metrics.xcLag};
snrs   = {spikes.clusters.metrics.snr};

I = ...
    ~cellfun(@isempty,  rpvs)  & ...
    ~cellfun(@isempty,  pAuc)  & ...
    ~cellfun(@isempty,ampRel)  & ...
    ~cellfun(@isempty, corrG)  & ...
    ~cellfun(@isempty, corrP)  & ...
    ~cellfun(@isempty,  subs)  & ...
    ~cellfun(@isempty, xcLag)  & ...
    ~cellfun(@isempty,  snrs);

I = find(I);

rpvs   = cell2mat(  rpvs(I));
pAuc   = cell2mat(  pAuc(I));
ampRel = cell2mat(ampRel(I));
corrG  = cell2mat( corrG(I));
corrP  = cell2mat( corrP(I));
subs   = cell2mat(  subs(I)); 
xcLag  = cell2mat( xcLag(I));
snrs   = cell2mat(  snrs(I));

X = [rpvs;pAuc;ampRel;corrG;corrP;subs;xcLag;snrs];
load('clusterClassifier.mat'); % load neural network
labels = net(X);

%% Initialize types to zero
nClusts = size(spikes.clusters.metrics,2);
for iClust = 1:nClusts; spikes.clusters.metrics(iClust).class = 0; end

%% Set cluster classes
nClusts = length(labels);
for iClust = 1:nClusts; spikes.clusters.metrics(I(iClust)).class = labels(iClust); end

end