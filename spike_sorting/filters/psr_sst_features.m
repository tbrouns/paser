function spikes = psr_sst_features(spikes,parameters)

spikes.features = [];
switch parameters.cluster.feature
    case 'pca'
        [~,features,~] = pca(psr_int16_to_single(spikes.waveforms(:,:),parameters));
        if (size(features,2) >= parameters.cluster.pca.dims)
            spikes.features = features(:,1:parameters.cluster.pca.dims)';
        end
    case 'wave'
        spikes = psr_sst_wavelet_features(spikes,parameters);
end
end