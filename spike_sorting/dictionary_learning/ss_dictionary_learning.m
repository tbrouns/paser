function assigns = ss_dictionary_learning(spikes,waveforms,fmm_p,align)

if (nargin < 3); fmm_p = spikes.params.fmm_p; end
if (nargin < 4); align = false; end

Sorter          = FMM(waveforms',spikes.params.fmm_k);
Sorter.align    = align;
Sorter.FMMparam = fmm_p;
Sorter.initialize;
Sorter.runVBfit;

assigns = getMAPassignment(Sorter);

end