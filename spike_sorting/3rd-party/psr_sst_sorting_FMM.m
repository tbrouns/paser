function spikes = psr_sst_sorting_FMM(spikes,parameters)

% Wrapper for dictionary learning spike sorting. 
% Focused mixture model (FMM) implementation.

addpath(genpath(parameters.path.fmm));

waveforms = psr_int16_to_single(spikes.waveforms(:,:),parameters); % Concatenate channels
waveforms = waveforms'; % Input should be [Npoints X Nspikes]

try 
    Sorter             = FMM(waveforms,parameters.sorting.fmm.k); 
    Sorter.align       = false;
    Sorter.FMMparam    = parameters.sorting.fmm.p;
    Sorter.maxClusters = parameters.sorting.fmm.kmax;
    Sorter.initialize;
    Sorter.runVBfit;
    spikes.assigns = getMAPassignment(Sorter);
catch ME
    str_1 = 'FMM ERROR:';
    str_2 = 'No spikes detected.';
    psr_show_warning({str_1,ME.message,str_2});
    spikes = [];
end

rmpath(genpath(parameters.path.fmm));

end