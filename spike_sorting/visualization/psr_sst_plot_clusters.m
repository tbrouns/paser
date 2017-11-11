function psr_sst_plot_clusters(spikes,clusters,savePath,filename)

if (nargin < 3); savePath = []; end % save in current working directory
if (nargin < 4); filename = []; 
else             filename = ['Clusters' filename(7:end)];
end

clustID = cell2mat({clusters.vars.id});
flags   = cell2mat({clusters.vars.flag});
clustID = clustID(flags);

nclusts = length(clustID);
nchan   = size(spikes.waveforms,3);

% Amplitude per channel

fig1 = figure; set(gcf,'position',get(0,'screensize'));
fig2 = figure; set(gcf,'position',get(0,'screensize'));
fig3 = figure; set(gcf,'position',get(0,'screensize'));
fig4 = figure; set(gcf,'position',get(0,'screensize'));

edges_min = -150:0;
edges_max =  0:150;

for iclust = 1:nclusts
    waves = spikes.waveforms(spikes.assigns == clustID(iclust),:,:);
    for ichan = 1:nchan
        amps_min_1 = min(waves(:,:,ichan),[],2);
        amps_max_1 = max(waves(:,:,ichan),[],2);
        
        figure(fig1);
        hold on
        subplot(nchan,1,ichan);
        histogram(amps_min_1,edges_min,'Normalization','probability')
        
        figure(fig2);
        hold on
        subplot(nchan,1,ichan);
        histogram(amps_max_1,edges_max,'Normalization','probability')
        
        for jchan = ichan:nchan
            
            kchan = ichan + (jchan - 1) * nchan;
            
            amps_min_2 = min(waves(:,:,ichan) + waves(:,:,jchan),[],2);
            amps_max_2 = max(waves(:,:,ichan) + waves(:,:,jchan),[],2);
            
            figure(fig3);
            hold on
            subplot(nchan,nchan,kchan);
            histogram(amps_min_2,edges_min,'Normalization','probability')
            title(['Channel ' num2str(ichan) ' vs ' num2str(jchan)]);
            
            figure(fig4);
            hold on
            subplot(nchan,nchan,kchan);
            histogram(amps_max_2,edges_max,'Normalization','probability')
            title(['Channel ' num2str(ichan) ' vs ' num2str(jchan)]);
        end
    end
end

figure(fig3); suptitle('Minimum cluster amplitudes');
figure(fig4); suptitle('Maximum cluster amplitudes');

% PCA

id    = cell2mat({clusters.vars.id});
flags = cell2mat({clusters.vars.flag});
id    = id(flags);

nclusts = length(id);

fig5 = figure; set(gcf,'position',get(0,'screensize'));

stds = zeros(nclusts,3);
avgs = zeros(nclusts,3);
for iclust = 1:nclusts
    which = find(spikes.assigns == id(iclust));
    x = spikes.waveforms(which,:) * spikes.info.pca.v(:,1);
    y = spikes.waveforms(which,:) * spikes.info.pca.v(:,2);
    z = spikes.waveforms(which,:) * spikes.info.pca.v(:,3);
    figure(fig5); hold on
    scatter3(x,y,z,'.');
    stds(iclust,:) = [ std(x), std(y), std(z)];
    avgs(iclust,:) = [mean(x),mean(y),mean(z)];
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