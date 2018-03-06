function spikes = psr_sst_sorting_CBP(spikes,data,parameters)

addpath(genpath(parameters.path.cbp));

% Create input structure

dataStruct.data       = data; clear data;
dataStruct.dt         = 1 / parameters.Fs;                
dataStruct.data       = psr_int16_to_single(dataStruct.data,parameters);
dataStruct.nchan      = size(dataStruct.data,1);
dataStruct.nsamples   = size(dataStruct.data,2);
dataStruct.polarity   = 'min'; % Change to parameter field
dataStruct.processing = [];

% Set parameters

params = load_default_parameters();
params.clustering.num_waveforms = parameters.sorting.cbp.nC;
params.general.plot_diagnostics = 0;
params.cbp.progress = false;
params.cbp.debug_mode = true;
filtdata = FilterData(dataStruct, params);

%% ----------------------------------------------------------------------------------
% Preprocessing Step 2: Estimate noise covariance and whiten

% Estimate and whiten the noise, assuming channel/time separability. This makes the
% L2-norm portion of the CBP objective into a sum of squares, simplifying the
% computation and improving computational efficiency.
%   - find "noise zones"
%   - compute noise auto-correlation function from these zones
%   - whiten each channel in time with the inverse matrix sqrt of the auto-correlation
%   - whiten across channels with inverse matrix sqrt of the covariance matrix.
% params.whitening includes:
%   - noise_threshold: noise zones must have cross-channel L2-norm less than this.
%   - min_zone_len: noise zones must have duration of at least this many samples.
%   - num_acf_lags: number of samples over which auto-correlation is estimated.

data_pp = WhitenNoise(filtdata, params);

%% ----------------------------------------------------------------------------------
% Preprocessing Step 3: Estimate initial spike waveforms

% Initialize spike waveforms, using clustering:
%  - collect data segments with L2-norm larger than params.clustering.spike_threshold
%  - align peaks of waveforms within these segments
%  - Perform PCA on these segments, select a subspace containing desired percent of variance
%  - Perform K-means clustering in this subspace
% params.clustering includes:
%  - num_waveforms : number of cells to be recovered
%  - spike_threshold : threshold used to pick spike-containing data segments (in stdevs)
%  - percent_variance : used to determine number of principal components to use for clustering

[centroids,~,~,~,~,~] = EstimateInitialWaveforms(data_pp, params);

init_waveforms = waveformMat2Cell(centroids, params.general.waveform_len, ...
                                  data_pp.nchan, params.clustering.num_waveforms);

% At this point, waveforms of all potential cells should be identified (NOTE:
% spike identification errors are irrelevant - only the WAVEFORMS matter).  If
% not, may need to adjust params.clustering.num_waveforms and re-run the clustering
% to identify more/fewer cells.  May also wish to adjust the
% params.general.waveform_len, increasing it if the waveforms (Fig 5) are being
% chopped off, or shortening it if there is a substantial boundary region of silence.
% If you do this, you should go back and re-run starting from the whitening step,
% since the waveform_len affects the identification of noise regions.

%% -----------------------------------------------------------------------------------
% CBP preprocessing

closeIfOpen(params.plotting.first_fig_num+[1,3]); % close Fourier and clustering figures

% To speed up computation, partition data into "snippets", which will be processed
% independently. Snippets have duration between min/max_snippet_len and are separated
% by "noise zones" in which the RMS of the waveforms does not surpass "threshold" for
% at least "min_separation_len" consecutive samples. Choose a conservative (low)
% threshold to avoid dropping spikes!

[snippets, ~, ~, snippet_centers, ~] = PartitionSignal(data_pp.data, params.partition);

%% ----------------------------------------------------------------------------------
% Use CBP to estimate spike times
 
[spike_times, ~, ~] = SpikesortCBP(snippets, snippet_centers, init_waveforms, params.cbp_outer, params.cbp);

nClusts = length(spike_times);
spiketimes = [];
assigns    = [];
for iClust = 1:nClusts
   times = spike_times{iClust};
   spiketimes = [spiketimes;times];    
   assigns    = [assigns;iClust*size(times)];
end
[spiketimes,I] = sort(spiketimes);
assigns = assigns(I);

spikes = psr_convert2spikes(spikes,data,spiketimes,assigns,parameters);

rmpath(genpath(parameters.path.cbp));

end