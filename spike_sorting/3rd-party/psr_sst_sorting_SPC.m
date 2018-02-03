function assigns = psr_sst_sorting_SPC(spikes,parameters)

nspikes = length(spikes.spiketimes);
mcs     = round(parameters.sorting.spc.mcs * nspikes);

inspk = psr_sst_wavelet_features(spikes,parameters); clear spikes;

%  * INPUT
%  *  D - Ndimensions-by-Npoints array of feature data
%  *  MCS - minimum cluster size (don't return clusters with <MCS points)
%  *
%  * OUTPUT
%  *  C - The clusters.  Cell array of #clusters vectors, where each vector 
%  *      contains a list of indexes of points in D which are in that cluster
%  *  P - Cluster parents.  Cluster i's parent cluster is C{P(i)}
% 
%  * SETTINGS
%  *  NMC - number of max-size clusters to consider
%  *  K - number of K-nearest neighbors to use

[C,P] = spc_mex(inspk,mcs);

% Convert to regular assigns

nclusts = length(C);

assigns = zeros(1,nspikes);
for iClust = 1:nclusts
    id = C{iClust} + 1;    
    assigns(id) = iClust;    
end
assigns = assigns + 1;

end