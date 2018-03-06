function spikes = psr_sst_sorting_FMM(spikes,parameters)

% Wrapper for dictionary learning spike sorting. 
% Focused mixture model (FMM) implementation.

addpath(genpath(parameters.path.fmm));

waveforms = psr_int16_to_single(spikes.waveforms(:,:),parameters); % Concatenate channels
waveforms = waveforms'; % Input should be [Npoints X Nspikes]

Sorter             = FMM(waveforms,parameters.sorting.fmm.k); 
Sorter.align       = false;
Sorter.FMMparam    = parameters.sorting.fmm.p;
Sorter.maxClusters = parameters.sorting.fmm.kmax;
Sorter.initialize;
Sorter.runVBfit;

spikes.assigns = getMAPassignment(Sorter);

rmpath(genpath(parameters.path.fmm));

end