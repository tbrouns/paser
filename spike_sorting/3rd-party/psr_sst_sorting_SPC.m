function spikes = psr_sst_sorting_SPC(spikes,parameters)

nSpikes  = length(spikes.spiketimes);
mcs      = round(spikes.info.dur * parameters.sorting.spc.mcs);
spikes   = psr_sst_features(spikes,parameters);
features = double(spikes.features);

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

[C,~] = spc_mex(features,mcs);

% Convert to regular assigns

nClusts = length(C);
assigns = zeros(1,nSpikes);
for iClust = 1:nClusts
    id = C{iClust} + 1;    
    assigns(id) = iClust;    
end
spikes.assigns = assigns + 1;

end