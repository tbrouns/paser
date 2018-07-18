function spikes = psr_sst_features(spikes,parameters)

spikes.features = [];
switch parameters.cluster.feature
    case 'pca'
        nChans = size(spikes.waveforms,3);
        for iChan = 1:nChans
            [~,features,~] = pca(psr_int16_to_single(spikes.waveforms(:,:,iChan),parameters));
            if (size(features,2) >= parameters.cluster.pca.dims)
                spikes.features = [spikes.features;features(:,1:parameters.cluster.pca.dims)'];
            end
        end
    case 'wave'
        spikes = psr_sst_wavelet_features(spikes,parameters);
end
end