function spikes = ums_clustering(spikes)

options.divisions = round(log2(1 / spikes.params.kmeans_clustersize));
options.divisions = max(min(options.divisions, spikes.params.divisions_max), spikes.params.divisions_min);

spikes = ss_align(spikes);
spikes = ss_kmeans(spikes,options);
spikes = ss_energy(spikes);

end