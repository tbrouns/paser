function assigns = ept_sst_sorting_SPC(spikes,parameters)

dims = parameters.sorting.spc.dims;
mcs  = parameters.cluster.size_min;
nspikes = length(spikes.spiketimes);

% Do PCA

PC   = pca(spikes.waveforms(:,:)');
if (size(PC,2) > dims)
   PC = PC(:,1:dims); % Take first D principle components
end

%  * INPUT
%  *  D - Ndimensions-by-Npoints array of feature data
%  *  MCS - minimum cluster size (don't return clusters with <MCS points)
%  *
%  * OUTPUT
%  *  C - The clusters.  Cell array of #clusters vectors, where each vector 
%  *      contains a list of indexes of points in D which are in that cluster
%  *  P - Cluster parents.  Cluster i's parent cluster is C{P(i)}

[C,P] = spc_mex(PC',mcs);

% Convert to regular assigns

nclusts = length(C);

assigns = zeros(1,nspikes);
for iClust = 1:nclusts
    id = C{iClust};    
    assigns(id) = iClust;    
end

end