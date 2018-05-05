function combis = psr_ft_combis(X,clustIDs)

nClusts  = length(clustIDs);
nCombis  = 0.5 * (nClusts^2 - nClusts);
combis   = cell(nCombis,2);
itr      = 1;
for iClust = 1:nClusts
    for jClust = iClust+1:nClusts
        combis(itr,:) = [{X.label{clustIDs(iClust)}},{X.label{clustIDs(jClust)}}];
        itr = itr + 1;
    end 
end

end