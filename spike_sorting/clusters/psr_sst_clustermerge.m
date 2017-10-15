function spikes = psr_sst_clustermerge(spikes,freq,parameters)

% Merge

if (~isfield(spikes,'assigns_prior'))
    spikes.assigns_prior = spikes.assigns; % Store prior-merge cluster assigns
end

if (strcmpi(parameters.sorting.method,'kst'))
    M = spikes.info.kst.simScore;
    M(~triu(true(size(M)),1)) = -realmax; % ignore lower triangular part
    IDs = find(M(:) >= parameters.cluster.merge_thresh);
else
    M = spikes.clusters.mahal;
    M(~triu(true(size(M)),1)) =  realmax; % ignore lower triangular part
    IDs = find(M(:) <= parameters.cluster.merge_thresh);
end

[I_row, I_col] = ind2sub(size(M),IDs);

clusterIDs = 1:size(M,1);
spikes = merge_clusters(spikes,clusterIDs(I_row),clusterIDs(I_col));
spikes = psr_sst_clusterfeatures(spikes,freq,parameters);

end

function spikes = merge_clusters(spikes,clustersX,clustersY)

nMerges = length(clustersX);
clustersXOld = clustersX;
clustersYOld = clustersY;

for iMerge = 1:nMerges
    
    % Cluster pair to merge
    
    clusterX = clustersX(iMerge);
    clusterY = clustersY(iMerge);
    
    % Find number of spikes for each cluster
    
    nSpikesX = sum((spikes.assigns == clusterX));
    nSpikesY = sum((spikes.assigns == clusterY));
        
    % Change cluster ID of smaller cluster to ID of larger cluster
    
    if (nSpikesX >= nSpikesY) 
        clustersX(clustersX == clusterY) = clusterX;
        clustersY(clustersY == clusterY) = clusterX;
    else
        clustersX(clustersX == clusterX) = clusterY;
        clustersY(clustersY == clusterX) = clusterY;
    end
    
end

[clustersXOld,IX] = unique(clustersXOld);
[clustersYOld,IY] = unique(clustersYOld);

clustersX = clustersX(IX);
clustersY = clustersY(IY);

for iMerge = 1:length(clustersX)
    spikes.assigns(spikes.assigns == clustersXOld(iMerge)) = clustersX(iMerge);
end

for iMerge = 1:length(clustersY)
    spikes.assigns(spikes.assigns == clustersYOld(iMerge)) = clustersY(iMerge);
end

end