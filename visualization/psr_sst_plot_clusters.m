function psr_sst_plot_clusters(spikes,parameters,savePath,filename)

if (nargin < 3); savePath = []; end % save in current working directory
if (nargin < 4); filename = []; 
else             filename = ['Clusters' filename(7:end)];
end

if (isfield(spikes,'delete')); spikes = psr_sst_filter_spikes(spikes,parameters,'delete'); end

% clusters = spikes.clusters;

% clustID = cell2mat({clusters.metrics.id});
% flags   = cell2mat({clusters.metrics.flag});
% clustID = clustID(flags);

% PCA

% id    = cell2mat({clusters.metrics.id});
% flags = cell2mat({clusters.metrics.flag});
% id    = id(flags);

id = [1 17];

nClusts = length(id);

fig5 = figure; set(gcf,'position',get(0,'screensize'));
figure(fig5); hold on

x = spikes.features(1,:);
y = spikes.features(2,:);
z = spikes.features(3,:);
scatter3(x,y,z,'.','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7]);
    
stds = zeros(nClusts,3);
avgs = zeros(nClusts,3);
for iClust = 1:nClusts
    spikeIDs = find(spikes.assigns == id(iClust));
    x = spikes.features(1,spikeIDs);
    y = spikes.features(2,spikeIDs);
    z = spikes.features(3,spikeIDs);
    scatter3(x,y,z,'.');
    stds(iClust,:) = [ std(x), std(y), std(z)];
    avgs(iClust,:) = [mean(x),mean(y),mean(z)];
end
stds = 3 * mean(stds);
avgs =     mean(avgs);
ux = avgs(1); uy = avgs(2); uz = avgs(3);
sx = stds(1); sy = stds(2); sz = stds(3);
axis([ux-sx ux+sx uy-sy uy+sy uz-sz uz+sz]);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

figure(fig1); export_fig([savePath filename '_P1']);
figure(fig2); export_fig([savePath filename '_P2']);
figure(fig3); export_fig([savePath filename '_P3']);
figure(fig4); export_fig([savePath filename '_P4']);
figure(fig5); view(3); export_fig([savePath filename '_P5']);
figure(fig5); view(2); export_fig([savePath filename '_P6']);

end