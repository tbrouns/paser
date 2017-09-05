function [spikes,clusters] = ss_clustermerge(spikes)

clusters = ss_clusterfeatures(spikes);

% only merge significant amplitude clusters

tf = ~cell2mat({clusters.vars.flag});
clusters.vars(tf)            = [];
clusters.mahal(tf,:)         = [];
clusters.mahal(:,tf)         = [];
clusters.bhattacharyya(tf,:) = [];
clusters.bhattacharyya(:,tf) = [];

% merge

spikes.assigns_prior = spikes.assigns;

% clustID  = cell2mat({clusters.vars.id});
% nclusts  = length(clustID);
% I_row = 1;
% 
% while (nclusts > 1 && ~isempty(I_row))
    
    M = clusters.mahal;
    
%     B = clusters.bhattacharyya;
    
    M(~triu(true(size(M)),1)) = realmax;
    
%     B(~triu(true(size(B)),1)) = realmax;

%     [mahal_min,I_M] = min(M(:));
%     [bhatt_min,I_B] = min(B(:));
    
    IDs = find(M(:) < spikes.params.thresh_mahal);
    [I_row, I_col] = ind2sub(size(M),IDs);
        
%     if (mahal_min < spikes.params.thresh_mahal && ~isempty(mahal_min))
%         [I_row, I_col] = ind2sub(size(M),I_M);
%         %     elseif (bhatt_min < spikes.params.thresh_bhattacharyya && ~isempty(bhatt_min))
%         %         [I_row, I_col] = ind2sub(size(B),I_B);
%     else
%         break;
%     end
        
    id = cell2mat({clusters.vars.id});
    spikes = merge_clusters(spikes,id(I_row),id(I_col));
        
    clusters = ss_clusterfeatures(spikes);
%     clustID  = cell2mat({clusters.vars.id});
%     nclusts  = length(clustID);
%     
% end

end

function spikes = merge_clusters(spikes,I1,I2)

N = length(I1);
I1_old = I1;
I2_old = I2;

for i = 1:N
    
    i1 = I1(i);
    i2 = I2(i);
    
    id1 = (spikes.assigns == i1);
    id2 = (spikes.assigns == i2);
    
    n1 = sum(id1);
    n2 = sum(id2);
    
    if (n1 >= n2);
        I1(I1 == i2) = i1;
        I2(I2 == i2) = i1;
    else
        I1(I1 == i1) = i2;
        I2(I2 == i1) = i2;
    end   
    
end

% id = (I1 == I2);
% I1(id) = [];
% I2(id) = [];
% N = length(I1);

for i = 1:N
        
    spikes.assigns(spikes.assigns == I1_old(i)) = I1(i);
    spikes.assigns(spikes.assigns == I2_old(i)) = I2(i);   
    
end


end