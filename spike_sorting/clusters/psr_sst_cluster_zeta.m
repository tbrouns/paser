function zetaDist = psr_sst_cluster_zeta(spikes,parameters)

% Find normalized zeta distance between clusters for merging [see Ref. 1, p.14]
% 
% References:
% [1] Yger, Pierre, et al. "Fast and accurate spike sorting in vitro and in
% vivo for up to thousands of electrodes." bioRxiv (2016): 067843.

% Filter spikes
spikes = psr_sst_filter_spikes(spikes,parameters,'delete');

nClust = max(spikes.assigns);
zetaDist = NaN(nClust,nClust);

for iClust = 1:nClust
    
    % Extract cluster ID
    
    for jClust = iClust+1:nClust
        
        spikeIDs_1 = ismember(spikes.assigns, iClust);
        spikeIDs_2 = ismember(spikes.assigns, jClust);
        
        if (~any(spikeIDs_1) || ~any(spikeIDs_2)); continue; end

        Y1 = spikes.features(:,spikeIDs_1);
        Y2 = spikes.features(:,spikeIDs_2);
                
        a1 = median(Y1,2);
        a2 = median(Y2,2);

        g = a1 - a2; 

        x1 = (g' * Y1);
        x2 = (g' * Y2);

        b1 = median(abs(x1 - median(x1))); % median absolute deviation
        b2 = median(abs(x2 - median(x2)));

        zetaDist(iClust,jClust) = norm(median(x1) - median(x2)) / (sqrt(b1^2 + b2^2));

    end
end

end