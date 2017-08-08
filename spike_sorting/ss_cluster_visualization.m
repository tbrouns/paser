function ss_cluster_visualization(spikes,clusters)

close all

clustID = cell2mat({clusters.vars.id});
flags   = cell2mat({clusters.vars.flag});
clustID = clustID(flags);

nclusts = length(clustID);
nchan   = size(spikes.waveforms,3);

% Amplitude per channel

fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;

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

            figure(fig4);
            hold on
            subplot(nchan,nchan,kchan);
            histogram(amps_max_2,edges_max,'Normalization','probability')
            
        end
    end
end

% PCA

id    = cell2mat({clusters.vars.id});
flags = cell2mat({clusters.vars.flag});
id    = id(flags);

nclusts = length(id);

figure;
hold on
for iclust = 1:nclusts
    which = find(spikes.assigns == id(iclust));
    x = spikes.waveforms(which,:) * spikes.info.pca.v(:,1);
    y = spikes.waveforms(which,:) * spikes.info.pca.v(:,2);
    z = spikes.waveforms(which,:) * spikes.info.pca.v(:,3);
    scatter3(x,y,z,'.');
end

end