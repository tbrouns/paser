function createSpikeHistogram(spikes,clusters)

tLength = spikes.info.detect.dur;
tBin    = 0.200; % sec
Fs      = spikes.params.Fs;

if (~isfield(spikes,'artifacts')); control = 1;
else                               control = 0;
end

clusterIDs = cell2mat({clusters.vars.id});
numclusts  = length(clusterIDs);

if (~control); artifactTimes = spikes.artifacts; end

X = 0:tBin:tLength+tBin;
x = X(1:end-1) + 0.5 * tBin;

for iClus = 1:numclusts

    id = (spikes.assigns == clusterIDs(iClus));
    spiketimes = spikes.spiketimes(id);

    N = binVector(spiketimes,X);
    
    figure;
    hold on
    bar(x,N);

    if (~control); scatter(artifactTimes,ones(size(artifactTimes))); end        
end

end

function N = binVector(times,X)

Y = discretize(times,X);
N = zeros(1,length(X)-1);
for i = 1:length(X)-1
   N(i) = sum(Y == i);
end

end