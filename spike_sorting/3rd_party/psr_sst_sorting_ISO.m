function assigns = psr_sst_sorting_ISO(spikes,parameters)

addpath(genpath(parameters.path.iso));

% Perform clustering using isotonic regression
nspikes = length(spikes.spiketimes);
dims    = parameters.sorting.wav.dims;
scales  = parameters.sorting.wav.scales;
mcs     = round(parameters.sorting.spc.mcs * nspikes);

opts = [];
opts.min_cluster_size = mcs;
opts.verbose = 0; % show plots (1: yes; 0: no) 
if (isfield(spikes,'assigns') && ~isempty(spikes.assigns))
    opts.initial_labels = spikes.assigns; 
end

waves = psr_single(spikes.waveforms(:,:),parameters); 
inspk = psr_wavelet_features(waves,dims,scales);
inspk = double(inspk');

[assigns,~] = isosplit5_mex(inspk,opts);

rmpath(genpath(parameters.path.iso));

end