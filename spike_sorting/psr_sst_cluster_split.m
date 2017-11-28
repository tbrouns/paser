function spikes = psr_sst_cluster_split(spikes,parameters)

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);

th(1,1,:) = spikes.info.thresh;
waves = psr_single(spikes.waveforms,parameters);
waves = waves ./ repmat(th, [size(waves,1) size(waves,2) 1]);

for iClust = 1:nClust
   
    clusterID = clusterIDs(iClust);
    which     = find(spikes.assigns == clusterID);
    ampClust  = max(waves(which,:),[],2);
    sd        = std(ampClust);
    bin       = sd * parameters.cluster.split_bin;
    edges     = min(ampClust):bin:max(ampClust);
    counts    = histcounts(ampClust,edges);
    counts    = gaussianSmoothing(counts,parameters.cluster.split_smooth);
    edges     = mean(diff(edges)) + edges(1:end-1);
    
    [pks,locs] = findpeaks(counts);
    [pks,I]    = sort(pks,'descend');
    locs       = locs(I);
    locs       = locs(1:2);
    pks        = min(pks(1:2));
    
    [~,I] = min(counts(locs(1):locs(2)));
    I = I + locs(1);
    ampDiff = pks - counts(I);
    
    if (ampDiff > parameters.cluster.split_thresh)
        id = ampClust > edges(I);
        id = which(id);
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