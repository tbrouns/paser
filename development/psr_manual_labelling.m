function labels = psr_manual_labelling(spikes,metadata,parameters,freq)

%------------- BEGIN CODE --------------

if (nargin < 4); freq = []; end

nClust = size(spikes.clusters.metrics,2);

% Convert and set necessary parameters

spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

spikes.info.kmeans.assigns = spikes.assigns;
numclusts                  = max(spikes.info.kmeans.assigns); % perhaps redundant? see nClust
cmap                       = jetm(numclusts);
spikes.info.kmeans.colors  = cmap(randperm(numclusts),:);

parameters = psr_parameters_display(parameters);
parameters.display.metrics       = false;
parameters.display.show_gaussfit = false;

%% Sort clusters based on similarity score

clusterIDs = [spikes.clusters.metrics.id];

% Get sim score
M = psr_sst_cluster_corr(spikes,parameters);
M = M(clusterIDs,:);
M = M(:,clusterIDs);

I  = triu(true(size(M)),1);
id = find(I);
[Ix,Iy] = ind2sub(size(I),id);

% Starting cluster
M = M(I);
[~,Imax] = max(M);
Ci = Ix(Imax);
Iz = [Ix,Iy];

% Initialize
clusterIDs = zeros(nClust,1);

for i = 1:nClust
    
    % Find next cluster index
    I = find(Iz(:,1) == Ci | Iz(:,2) == Ci);
    [~,id] = max(M(I));
    Imax = I(id);
    Iz_row = Iz(Imax,:);
    Cj = Iz_row(Iz_row~=Ci);
    clusterIDs(i) = Ci;
    
    % Remove current cluster index from arrays
    id = (Iz(:,1) == Ci | Iz(:,2) == Ci);
    Iz(id,:) = [];
    M (id)   = [];
    
    Ci = Cj;
end

ids = [spikes.clusters.metrics.id];
clusterIDs = ids(clusterIDs);

run(parameters.general.configPath); % TEMP

%% Ignore low spike count clusters
nspikes = zeros(nClust,1);
for iClust = 1:nClust
    nspikes(iClust) = sum(spikes.assigns == clusterIDs(iClust));
end
I = nspikes > parameters.cluster.quality.min_spikes;
clusterIDs = clusterIDs(I);
nClust = length(clusterIDs);

%% Manual labelling

fig = figure; set(gcf,'position',get(0,'screensize'));

labels = zeros(nClust,1);
iClust = 1;
while iClust <= nClust
    
    clusterID = clusterIDs(iClust);
    
    figure(fig); clf;
    subaxis(2,3,1:2,'Margin',0.05,'Padding',0); psr_sst_plot_waveforms(spikes,clusterID,parameters);
    subaxis(2,3,  3,'Margin',0.05,'Padding',0); psr_sst_plot_amp (spikes,clusterID,parameters);
    subaxis(2,3,4:5,'Margin',0.05,'Padding',0); psr_sst_plot_stability(spikes,clusterID,freq,metadata,parameters);
    subaxis(2,3,  6,'Margin',0.05,'Padding',0); psr_sst_plot_isi      (spikes,clusterID,parameters);
    
    % Wait for key press
    but = waitForButton();
    
    % Process key press
    
    id = find(ids == clusterID);
    if     isequal(but,49); labels(id) = 1; % 1
    elseif isequal(but,50); labels(id) = 2; % 2
    elseif isequal(but,51); labels(id) = 3; % 3
    elseif isequal(but,52); labels(id) = 4; % 4
    elseif isequal(but,53); labels(id) = 5; % 5
    elseif isequal(but,29) % right arrow
        % skip
    elseif isequal(but,28) % left arrow
        iClust = iClust - 2; % go back to previous clusters
    else % stay with current cluster
        iClust = iClust - 1;
    end
    
    iClust = iClust + 1;
    if iClust < 1; iClust = 1; 
    elseif iClust > nClust
        suptitle('Press key to continue');
        but = waitForButton();
        if isequal(but,28) % left arrow
            iClust = nClust;
        end
    end
    
end

close(fig);

end

function but = waitForButton()

m = 0;
while ~m
    m   = waitforbuttonpress;
    but = double(get(gcf,'CurrentCharacter'));
end

end

%%

%------------- END OF CODE --------------