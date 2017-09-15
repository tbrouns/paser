function spikes = ss_spikefilter_pca(spikes)

% Remove outliers based on PC distance

for iclust = 1:nclusts
    
    which = find(spikes.assigns == clustID(iclust));
    x = spikes.waveforms(which,:) * spikes.info.pca.v(:,1);
    y = spikes.waveforms(which,:) * spikes.info.pca.v(:,2);
    z = spikes.waveforms(which,:) * spikes.info.pca.v(:,3);

    M  = [x,y,z];
    D  = mahal(M,M);
    sd = std(D);
    id = D > spikes.params.outlier.std * sd;
    id = which(id);
    
    spikes = ss_spike_removal(spikes,id);
    
end

end