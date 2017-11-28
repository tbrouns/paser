function assigns = psr_sst_sorting_FMM(spikes,parameters)

% Wrapper for dictionary learning spike sorting. 
% Focused mixture model (FMM) implementation.

waveforms = spikes.waveforms(:,:); % Concatenate channels
waveforms = waveforms'; % Input should be [Nsamples X Nspikes]

Sorter             = FMM(waveforms,parameters.sorting.fmm.k); 
Sorter.align       = false;
Sorter.FMMparam    = parameters.sorting.fmm.p;
Sorter.maxClusters = parameters.sorting.fmm.kmax;
Sorter.initialize;
Sorter.runVBfit;

assigns = getMAPassignment(Sorter);

end