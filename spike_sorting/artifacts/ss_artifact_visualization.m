function ss_artifact_visualization(spikes)

close all

waveforms = spikes.artifact_waveforms;
waveforms = reshape(waveforms,[],size(waveforms,2) * size(waveforms,3));

assigns_all = ss_dictionary_learning(spikes,waveforms);

nclusters = length(unique(assigns_all));

fig1 = figure;
fig2 = figure;
t_step = 0.001; % s
edges = 0:t_step:120;

for iclust = 1:nclusters
    clf
    id = assigns_all == iclust;
    waves = waveforms(id,:);
    waveforms_mean = mean(waveforms);
    t = 1:size(waveforms_mean,2);    
    
    figure(fig1);
    cmap    = hot(64);
    [n,x,y] = histxt(waves);
    h       = imagesc(x,y,n);
    colormap(cmap);
    set(gca,'Color', cmap(1,:) );
    xlabel('Sample #');
    ylabel('Voltage');
        
    figure(fig2);
    artifactTimes = spikes.artifacts(id);
    histogram(artifactTimes,edges)
    xlabel('Time [s]');
    ylabel('Spike count');
    title(['Bin size = ' num2str(t_step) ' sec']);
    
end

end