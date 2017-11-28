function spikes = ss_energy(spikes)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_energy - Interface energy based cluster similarity computation.
%
% Usage:
%      spikes = ss_energy(spikes)
%
% Description:
%     SPIKES = SS_ENERGY(SPIKES) adds an interface energy matrix to a
%     spike-sorting object in SPIKES.INFO.INTERFACE_ENERGY.
%
%     The energy similarity matrix is calculated by applying an exponential
%     decay to all pairwise euclidean distances between waveforms from two
%     clusters (or within a single cluster for intra-cluster energy) and
%     summing these distances.
%
%     The calculation ignores the energy due to the zero distance between
%     points and themselves; this removes a dependence of the density on
%     the absolute size of the cluster.  As a result, singleton clusters
%     do not have a well-defined energy and will cause an error.
%
%     When each entry is normalized by the number of distinct contributing
%     pairs (Na*Nb for off diagonal entries and Na*(Na-1)/2 on the diagonal),
%     it approximates the fraction of pairs in a given cluster whose distance
%     is not much greater than the length constant of the exponential and thus
%     provides an estimate of local density.  This function does not, however,
%     normalize SPIKES.INFO.INTERFACE_ENERGY, since the normalized form is
%     inconvenient during cluster aggregation.  The normalization can readily
%     be done, however, with
%          normalize = ((numpts * numpts') - diag(numpts));
%          normalize = normalize - diag(0.5*diag(normalize));
%          normalized_energy = interface_energy ./ normalize;
%     where 'numpts' is a vector of cluster sizes.
%
%     The unnormalized energy matrix can be updated during aggregation without
%     the need to recompute it from scratch.  The intra-cluster energy E(AB,AB)
%     of a cluster AB formed by aggregating clusters A and B is given by
%              E(AB,AB) = E(A,A) + E(B,B) + E(A,B)
%     and the inter-cluster energy between any cluster C and an aggregate AB is
%                 E(AB,C) = E(A,C) + E(B,C)
%
%      Note that this function also renumbers the miniclusters based on
%      connection strength.
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%
% Major edits made by Terence Brouns (2017).

%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') || (size(spikes.waveforms, 1) < 1))
    disp('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
    return;
elseif (~isfield(spikes.info, 'kmeans'))
    disp('SS:overcluster_not_computed', 'The data must be overclustered before computing energy');
    return
end

nClusts = length(unique(spikes.info.kmeans.assigns));
d       = diag(spikes.info.pca.s);
r       = find(cumsum(d)/sum(d) > 0.95, 1);
waves   = (spikes.waveforms(:,:) * spikes.info.pca.v(:,1:r))';

%%%%% PREPARE SOME INFORMATION

waveclust = cell(nClusts,1); % collect spikes for each cluster
numpts    = zeros(nClusts,1);
for clust = 1:nClusts
    waveclust{clust} = waves(:,spikes.info.kmeans.assigns == clust);
    numpts(clust)    = size(waveclust{clust},2);
end
clear waves

%%%%% HEURISTIC DISTANCE SCALE that seems to work.  The calculation is not too sensitive to this parameter.
scale = sqrt(sum(diag(spikes.info.kmeans.W))) / 10;

%%%%% PREPARE TO LOOP

k       = 0;
sumpts  = sum(tril(numpts * ones(1,length(numpts))));
total   = sum(sumpts);
progress_bar(0,1,'Computing Interaction Energies . . .')
interface_energy  = zeros(nClusts);

%%%%% PAIRWISE DISTANCES LOOP

for iClust = 1:nClusts
    X  = waveclust{iClust};
    for jClust = iClust:nClusts   % clust2 starts at clust1 so we get intra- too
        % Compute pairwise Euclidean distance matrix
        % use formula that (x-y)^2 = x^2 + y^2 - 2xy
        % Note that this formula can cause small negative values due to
        % round-off, so we take absolute value
        Y = waveclust{jClust};
        dists = abs(bsxfun(@plus,dot(X,X,1),dot(Y,Y,1)')-(2*Y'*X));
        interface_energy(iClust,jClust) = sum(exp(-realsqrt(dists(:))/scale));
    end
    k = k + sumpts(iClust);
    progress_bar(k/total,[]);
end

%%%%% CORRECTION TERMS
% The energy matrix so far includes a contribution in the intra-cluster
% energies that is not found in the inter-cluster energies; namely, the
% computation of   sum_(all x) sum_(all y) e^(-dist/scale)   for
% intra-cluster energy includes cases where x == y (so dist == 0).
interface_energy = interface_energy - diag(numpts);     % So subtract this out.

% Also, we've double counted pairs in the intra-energy case, since dist(a,b)
% and dist(b,a) are not treated as distinct;
interface_energy = interface_energy - diag(0.5*diag(interface_energy));

%%%%% FINISH UP
spikes.info.interface_energy = interface_energy;

% Now re-number all clusters based on connection strength similarities
assignments      = double(spikes.info.kmeans.assigns);
interface_energy = spikes.info.interface_energy;
nClusts        = max(assignments);
numpts           = full(sparse(assignments, 1, 1, nClusts, 1));
normalize        = ((numpts * numpts') - diag(numpts));           % Off diag: Na*Nb, On diag: Na^2-Na ...
normalize        = normalize - diag(0.5 * diag(normalize));       % ... and divide diagonal by 2
norm_energy      = interface_energy ./ normalize;
self             = repmat(diag(norm_energy), [1,nClusts]);
connect_strength = 2 .* norm_energy ./ (self + self');
connect_strength = connect_strength .* (1-eye(nClusts));        % diag entries <- 0, so we won't agg clusters with themselves

% initialize
[~,pos(1)] = max(max(connect_strength));
for j = 2:nClusts
    [val1,pos1] = max(connect_strength(:,pos(j-1)));
    [val2,pos2] = max(connect_strength(pos(j-1),:));
    if val1 >= val2; pos(j) = pos1;
    else             pos(j) = pos2;
    end
    connect_strength(:,pos(j-1)) = -inf;
    connect_strength(pos(j-1),:) = -inf;
end

% update kmeans assignments
assigns = spikes.info.kmeans.assigns;
ie      = spikes.info.interface_energy;
c       = spikes.info.kmeans.centroids;

for j = 1:nClusts
    assigns(spikes.info.kmeans.assigns == pos(j)) = j;
    for k = 1:nClusts
        if k > j
            if pos(k) > pos(j); ie(j,k) = spikes.info.interface_energy(pos(j),pos(k));
            else                ie(j,k) = spikes.info.interface_energy(pos(k),pos(j));
            end
        end
    end
    c(j,:) = spikes.info.kmeans.centroids(pos(j),:);
end

spikes.info.kmeans.assigns   = assigns;
spikes.info.interface_energy = ie;
spikes.info.kmeans.centroids = c;

