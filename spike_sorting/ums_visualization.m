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
    
    if length(find(spikes.assigns==clus)) > 10 % visualize the cluster only if there are more than 10 spikes within
        clf
        subplot(2,3,1); plot_waveforms(spikes,clus);
        subplot(2,3,2); plot_residuals(spikes,clus);
        subplot(2,3,4); plot_detection_criterion(spikes,clus);
        subplot(2,3,5); plot_isi(spikes,clus);
        subplot(2,3,6); plot_stability(spikes,clus);
        subplot(2,3,3); plot_distances(spikes,clus);
        xlim([0 200]);
        xlabel('Z-score'); ylabel('Count');
        export_fig([...
            'Spikes_'   metadata.subject ...
            '_'         metadata.session ...
            '_T'        num2str(metadata.tetrode,           '%02d')   ...
            '_AG'       num2str(spikes.params.agg_cutoff,   '%10.2e') ...
            '_ST'       num2str(spikes.params.thresh,       '%04.1f') ...
            '_Cluster'  num2str(i,                          '%03d')]);
        xcorr = cross_correlation(spikes,clus);    
        disp('Cross-correlation:');
        disp(xcorr);
    end
end
end

function xc = cross_correlation(spikes, show)

% calculate standard deviation
clus         = get_spike_indices(spikes, show );
memberwaves  = spikes.waveforms(clus,:);
num_channels = size(spikes.waveforms,3);
num_samples  = size(spikes.waveforms,2);
s            =  std(memberwaves);

% calculate cross-correlation

xc = zeros(num_channels); % cross-correlation matrix
sd = reshape(s,num_samples,num_channels); % residuals for each channel separately

for ichan = 1:num_channels
    for jchan = ichan:num_channels
        xc(ichan,jchan) = xcorr(sd(:,ichan),sd(:,jchan),0); % xcorr with zero lag
    end
end

waves = reshape(mean(memberwaves),num_samples,num_channels);  
id = double(max(abs(waves)) > abs(spikes.info.detect.thresh)); % check which channels are active
id = id' * id; 
xc = xc .* id; % ignore channels that do not display spike

end