function spikes = psr_sst_sorting_UMS(spikes,parameters)

% Data conversion

spikes.params.Fs                 = parameters.Fs;
spikes.params.detect.thresh      = parameters.spikes.thresh;
spikes.params.detect.window_size = parameters.spikes.window_size; 

spikes.info.detect.thresh = spikes.info.thresh;
spikes.info.detect.dur    = spikes.info.dur;
spikes.info.detect.stds   = spikes.info.stds;

% SVD the data matrix
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);
[pca.u,pca.s,pca.v] = svd(detrend(spikes.waveforms(:,:),'constant'), 0); 
spikes.info.pca = pca;

% Mean statistics
spikes.info.detect.stds   = mean(spikes.info.detect.stds);
spikes.info.detect.thresh = mean(spikes.info.detect.thresh);

spikes.params.agg_cutoff         = parameters.sorting.ums.agg_cutoff; 
spikes.params.kmeans_clustersize = parameters.sorting.ums.kmeans_size;

% Run through pipeline
spikes = ss_kmeans   (spikes);
spikes = ss_energy   (spikes);
spikes = ss_aggregate(spikes);

% Convert back to int16
precision = 10^parameters.general.precision;
spikes.waveforms = int16(precision * spikes.waveforms);
 

end