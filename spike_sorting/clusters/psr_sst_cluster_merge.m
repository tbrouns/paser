function spikes = psr_sst_cluster_merge(spikes,parameters)

% PSR_SST_CLUSTER_MERGE - Merge spike clusters based on waveform similarity.
% This function merges spike clusters obtained by spike sorting if their
% calculated similarity score is above a specified threshold. The metric
% used for the similarity score depends on the chosen method.
%
% Syntax:  spikes = psr_sst_cluster_merge(spikes,freq,metadata,parameters)
%
% Inputs:
%    spikes     - See README
%    freq       - See README
%    metadata   - See README
%    parameters - See README
%
% Outputs:
%    spikes - The field "spikes.assigns" contains the new cluster ID for
%             each spike after merging, while "spikes.assigns_prior" retains
%             the old cluster IDs before merging.
%
% See also: PSR_SST_CLUSTER_FEATURES

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

% Store prior-merge cluster assigns

if (~isfield(spikes,'assigns_prior')); spikes.assigns_prior = spikes.assigns; end
if (~isfield(spikes,'features')); disp('Spike features missing. Exiting function.'); return; end

spikes.assigns = spikes.assigns_prior;
nClust = max(spikes.assigns);

if (nClust > 0)
    
    clusterIDs = 1:nClust;
    
    while true
        
        [I_row, I_col, score] = findClustersToMerge(spikes,parameters);
        
        if (~isempty(I_row))
            
            %% Do the merging
            
            % Cluster pair to merge
            
            clusterX  = clusterIDs(I_row);
            clusterY  = clusterIDs(I_col);
            
            % Find number of spikes for each cluster
            
            nSpikesX = sum(spikes.assigns == clusterX);
            nSpikesY = sum(spikes.assigns == clusterY);
            
            % Change cluster ID of smaller cluster to ID of larger cluster
            
            if (nSpikesX >= nSpikesY); spikes.assigns(spikes.assigns == clusterY) = clusterX;
            else,                      spikes.assigns(spikes.assigns == clusterX) = clusterY;
            end
            
        else
            break;
        end
    end
    
    switch parameters.cluster.merge.type % Save data
        case 'zeta'; spikes.clusters.zeta = score;
        case 'corr'; spikes.clusters.corr = score;
        case 'bhat'; spikes.clusters.bhat = score;
    end
    
end

end

function [I_row, I_col, score] = findClustersToMerge(spikes,parameters)

score = [];
I_row = [];
I_col = [];

% Find which clusters to merge

switch parameters.cluster.merge.type
    case 'zeta'; score =  psr_sst_cluster_zeta(spikes,parameters); tf = score(:) <=  parameters.cluster.merge.zeta_thresh;
    case 'corr'; score = -psr_sst_cluster_corr(spikes,parameters); tf = score(:) <= -parameters.cluster.merge.corr_thresh;
    case 'bhat'; score =  psr_sst_cluster_bhat(spikes,parameters); tf = score(:) <=  parameters.cluster.merge.bhat_thresh;
end

[~,mergeID] = min(score(:));
if (tf(mergeID)); [I_row, I_col] = ind2sub(size(score),mergeID); end

end


%------------- END OF CODE --------------