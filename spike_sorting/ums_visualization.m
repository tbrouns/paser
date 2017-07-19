function ums_visualization(spikes,clusters,metadata,type)

if (nargin < 4); type = 'single'; end

close all

% Only plot single units

if (strcmp(type,'single'))
    tf = ~strcmp({clusters.vars(:).unit},{'single'});
    clusters.vars(tf) = [];
end

nclus = size(clusters.vars,2);

% Plot

fig1 = figure; set(gcf,'position',get(0,'screensize'));

for i = 1 : nclus
    
    icluster = clusters.vars(i).id;
    
    if length(find(spikes.assigns == icluster)) > 10 % visualize the cluster only if there are more than 10 spikes within
        figure(fig1);
        clf
        subplot(2,3,1); plot_waveforms(spikes,icluster);
        subplot(2,3,2); plot_residuals(spikes,icluster);
        subplot(2,3,3); plot_distances(spikes,icluster);
        subplot(2,3,4); plot_detection_criterion(spikes,icluster);
        subplot(2,3,5); plot_isi(spikes,icluster);
        subplot(2,3,6); plot_stability(spikes,icluster);
        export_fig([...
            'Spikes_'   metadata.subject ...
            '_'         metadata.session ...
            '_T'        num2str(metadata.tetrode,         '%02d')   ...
            '_'         spikes.params.cluster_method                ...
            '_AG'       num2str(spikes.params.agg_cutoff, '%10.2e') ...
            '_ST'       num2str(spikes.params.thresh,     '%04.1f') ...
            '_Cluster'  num2str(icluster,                 '%03d')]);
    end
end

end