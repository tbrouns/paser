function assigns = ept_sst_sorting_UMS(spikes,parameters)

spikes.params.detect.method      = parameters.spikes.method;
spikes.params.detect.thresh      = parameters.spikes.thresh;
spikes.params.detect.window_size = parameters.spikes.window_size; 
spikes.params.detect.shadow      = parameters.spikes.shadow;      
spikes.params.detect.cross_time  = parameters.spikes.cross_time;  
spikes.params.detect.max_jitter  = parameters.spikes.max_jitter; 

samples_before = round(spikes.params.Fs * spikes.params.detect.cross_time  / 1000);
jitter_range   = samples_before - 1 + (1:round(spikes.params.detect.max_jitter * spikes.params.Fs / 1000));

% identify which channel the event occurred on
divisor = repmat(spikes.info.detect.thresh, [size(spikes.waveforms,1) 1]);
[~, spikes.info.detect.event_channel] = max(squeeze(min(spikes.waveforms(:,jitter_range,:), [], 2))./divisor, [], 2);
spikes.info.detect.event_channel = single(spikes.info.detect.event_channel);

% SVD the data matrix
[pca.u,pca.s,pca.v] = svd(detrend(spikes.waveforms(:,:),'constant'), 0); 
spikes.info.pca = pca;

% Mean statistics
spikes.info.detect.stds   = mean(spikes.info.detect.stds);
spikes.info.detect.thresh = mean(spikes.info.detect.thresh);

% Convert parameters to UMS format

% UltraMegaSort2000
spikes.params.agg_cutoff         = parameters.sorting.ums.agg_cutoff; 
spikes.params.kmeans_clustersize = parameters.sorting.ums.kmeans_size;

% Run through pipeline
spikes = ss_kmeans   (spikes);
spikes = ss_energy   (spikes);
spikes = ss_aggregate(spikes);
assigns = spikes.assigns;
            
end