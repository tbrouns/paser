function spikes = psr_sst_cluster_split(spikes,parameters)

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);

th(1,1,:) = spikes.info.thresh;
waves = psr_single(spikes.waveforms,parameters);
waves = waves ./ repmat(th, [size(waves,1) size(waves,2) 1]);

for iClust = 1:nClust
   
    clusterID = clusterIDs(iClust);
    spikeIDs  = find(spikes.assigns == clusterID);
    ampClust  = max(waves(spikeIDs,:),[],2);
    sd        = std(ampClust);
    bin       = sd * parameters.cluster.split.bin;
    edges     = min(ampClust):bin:max(ampClust);
    counts    = histcounts(ampClust,edges); % amplitude distribution
    counts    = gaussianSmoothing(counts,parameters.cluster.split.smooth);
    edges     = mean(diff(edges)) + edges(1:end-1);
    
    [pks,locs] = findpeaks(counts);
    [pks,I]    = sort(pks,'descend');
    locs       = locs(I);
    locs       = locs(1:2); % Largest two peaks
    pks        = min(pks(1:2)); % Second largest peaks
    
    [~,I] = min(counts(locs(1):locs(2))); % Minimum in-between peaks
    I = I + locs(1);
    ampDiff = pks - counts(I); % Prominence of second largest peak 
    
    if (ampDiff > parameters.cluster.split.thresh)
        id = ampClust > edges(I);
        id = spikeIDs(id);
        maxId = max(unique(spikes.assigns));
        spikes.assigns(id) = maxId + 1;
    end
end

end

function x = gaussianSmoothing(x,n)

g = gausswin(n); % create Gaussian smoothing window
g = g / sum(g); % normalize
x = conv(x, g, 'same');

end