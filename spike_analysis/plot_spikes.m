function plot_spikes(spikes,clusters,type)

if strcmp(type,'single')
    tf = ~strcmp({clusters.vars(:).unit},{'single'});
    clusters.vars(tf) = [];
end

clustIds = cell2mat({clusters.vars.id});
numclusts = length(clustIds);

fig1 = figure;
fig2 = figure;

bins = 0:0.001:0.5;

for iclust = 1:numclusts
    figure(fig1);
    hold on
    spike_times = spikes.spiketimes(spikes.assigns == clustIds(iclust));
%     scatter(spike_times,clustIds(iclust)*ones(size(spike_times)),'.');
    
    % Draw spikes (T is 1xN matrix with N spikes at times t)
    height = 0.1;
    plot([spike_times;spike_times],[ones(size(spike_times))*(clustIds(iclust)-height).';ones(size(spike_times))*(clustIds(iclust)+height).'],'k');

    figure(fig2);
    subplot(numclusts,1,iclust);
    isi = diff(spike_times);
    h = histogram(isi,bins);
    h.FaceColor = [0 0 1.0];
    h.FaceAlpha = 0.9;
    h.EdgeColor = 'none';

end

figure(fig1);
ylim([min(clustIds) - 1, max(clustIds) + 1]);
xlabel('Time [s]');

end