function spikes = psr_sst_sorting_ISO(spikes,parameters)

% Perform clustering using isotonic regression
mcs = round(spikes.info.dur * parameters.sorting.iso.mcs);
 
opts = [];
opts.min_cluster_size = mcs;
opts.verbose = 0; % show plots (1: yes; 0: no) 
if (isfield(spikes,'assigns') && ~isempty(spikes.assigns))
    opts.initial_labels = spikes.assigns; 
end

spikes = psr_sst_features(spikes,parameters);

% [assigns,~] = isosplit5(inspk,opts);
features = double(spikes.features);
spikes.assigns = isosplit5_mex(features);

end