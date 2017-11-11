function [spikes,I] = psr_sst_cluster_sort(spikes)

clustIDs = unique(spikes.assigns);
nclusts  = length(clustIDs);
nspikes  = zeros(1,nclusts);
for iClust = 1:nclusts
    nspikes(iClust) = sum(spikes.assigns == clustIDs(iClust));
end

[~,I] = sort(nspikes,'descend');
r = 1:length(nspikes);
r(I) = r; % rank

assignsNew = spikes.assigns;
for iClust = 1:nclusts
   assignsNew(spikes.assigns == clustIDs(iClust)) = r(iClust);
end

spikes.assigns = assignsNew;

end