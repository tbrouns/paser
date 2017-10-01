function assigns = ept_sst_sorting_FMM(spikes,waveforms,parameters)

% Wrapper for dictionary learning spike sorting. Focused mixture model
% (FMM) implementation.

Sorter             = FMM(waveforms',spikes.params.sorting.fmm.k);
Sorter.align       = align;
Sorter.FMMparam    = parameters.sorting.fmm.p;
Sorter.maxClusters = parameters.sorting.fmm.kmax;
Sorter.initialize;
Sorter.runVBfit;

assigns = getMAPassignment(Sorter);

end