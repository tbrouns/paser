function ums_visualization(spikes,metadata)

close all

id          = unique(spikes.assigns);
nclusters   = length(id);
nwaves      = zeros(nclusters,2);
for icluster = 1:nclusters
    nwaves(icluster,1) = sum(spikes.assigns == id(icluster));
    nwaves(icluster,2) = id(icluster);
end

nwaves  = sortrows(nwaves,-1);
nclus   = size(nwaves,1);

figure;
set(gcf,'position',get(0,'screensize'));
   
for i = 1 : nclus
    clus = nwaves(i,2);
    
    if length(find(spikes.assigns==clus))>10 % visualize the cluster only if there are more than 10 spikes within
        clf
        subplot(2,3,1); plot_waveforms(spikes,clus);
        subplot(2,3,2); plot_residuals(spikes,clus);
        subplot(2,3,4); plot_detection_criterion(spikes,clus);
        subplot(2,3,5); plot_isi(spikes,clus);
        subplot(2,3,6); plot_stability(spikes,clus);
        subplot(2,3,3); plot_distances(spikes,clus);
        xlim([0 200]);
        xlabel('Z-score'); ylabel('Count');
        export_fig(['Clusters_' metadata.subject '_' metadata.session '_T' num2str(metadata.tetrode) '_AG' num2str(spikes.params.agg_cutoff) '_ST' num2str(spikes.params.thresh) '_Cluster' num2str(i)]);
    end
end
end