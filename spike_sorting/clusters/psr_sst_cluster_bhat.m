function bhattaDist = psr_sst_cluster_bhat(spikes,parameters)

% PSR_SST_CLUSTER_BHAT - Bhattacharyya distance between each pair of clusters
%
% Syntax:  bhattaDist = psr_sst_cluster_bhat(spikes,parameters)
%
% Inputs:
%    spikes     - See README
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    bhattaDist - Matrix of Bhattacharyya distances between every pair of
%                 clusters, with shape:
%                 [Number of clusters x Number of clusters]
%
% See also: PSR_SST_CLUSTER_MERGE

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Filter spikes
spikes = psr_sst_filter_spikes(spikes,parameters,'delete');

nClust = max(spikes.assigns);
bhattaDist = NaN(nClust,nClust); % Bhattacharyya distance

for iClust = 1:nClust
        
    for jClust = iClust+1:nClust
        
        spikeIDs_1 = ismember(spikes.assigns, iClust);
        spikeIDs_2 = ismember(spikes.assigns, jClust);
        
        if (~any(spikeIDs_1) || ~any(spikeIDs_2)); continue; end

        f1 = spikes.features(:,spikeIDs_1)';
        f2 = spikes.features(:,spikeIDs_2)';
                
        if (size(f1,1) > size(f1,2) && size(f2,1) > size(f2,2)) % Need to have more observations than dimensions
            try    bhattaDist(iClust,jClust) = bhattacharyya(f1,f2);
            catch; continue;
            end
        end
    end
end

end