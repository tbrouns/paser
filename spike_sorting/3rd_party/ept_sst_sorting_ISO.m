function assigns = ept_sst_sorting_ISO(spikes,parameters)

% Perform clustering using isotonic regression
dims = parameters.sorting.iso.dims;

opts = [];
opts.min_cluster_size = parameters.cluster.size_min;
opts.verbose = 0; % show plots (1: yes; 0: no) 

if (isfield(spikes,'assigns')); opts.initial_labels = spikes.assigns; end

% Do PCA

PC = pca(spikes.waveforms(:,:)');
if (size(PC,2) > dims)
   PC = PC(:,1:dims); % Take first D principle components
end

data = double(PC');

[assigns,~] = isosplit5_mex(data,opts);

end