function spikes = psr_manual_labelling(spikes,metadata,parameters,freq)

%------------- BEGIN CODE --------------

if (nargin < 4); freq = []; end

nClust = size(spikes.clusters.vars,2);

% Convert and set necessary parameters

if (isa(spikes.waveforms,'int16'))
    spikes.waveforms = psr_single(spikes.waveforms,parameters);
end

spikes = psr_sst_display_parameters(spikes);
spikes.unwrapped_times          = spikes.spiketimes;
spikes.params.detect.ref_period = parameters.spikes.ref_period;
spikes.params.detect.shadow     = 0.5 * parameters.spikes.window_size;
spikes.info.kmeans.assigns      = spikes.assigns;
numclusts                       = max(spikes.info.kmeans.assigns); % perhaps redundant? see nClust
cmap                            = jetm(numclusts);
spikes.info.kmeans.colors       = cmap(randperm(numclusts),:);


spikes.params.display.metrics = false;
spikes.params.display.show_gaussfit = false;

spikes = psr_sst_filter_amp (spikes,parameters,'array');
spikes = psr_sst_spike_removal(spikes,find(spikes.removed),'delete');

%% Sort clusters based on similarity score

clusterIDs = [spikes.clusters.vars.id]; 

% Get sim score
M = spikes.info.kst.simScore;
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

ids = [spikes.clusters.vars.id]; 
clusterIDs = ids(clusterIDs);

%% Manual labelling

fig = figure; set(gcf,'position',get(0,'screensize'));

iClust = 1;
while iClust <= nClust
    clusterID = clusterIDs(iClust);
    if length(find(spikes.assigns == clusterID)) > parameters.cluster.min_spikes
        figure(fig); clf;
        subaxis(2,3,1:2,'Margin',0.05,'Padding',0); psr_sst_plot_waveforms  (spikes,clusterID,parameters);
        subaxis(2,3,  3,'Margin',0.05,'Padding',0); plot_detection_criterion(spikes,clusterID,parameters);
        subaxis(2,3,4:5,'Margin',0.05,'Padding',0); psr_sst_plot_stability  (spikes,clusterID,freq,metadata,parameters);
        subaxis(2,3,  6,'Margin',0.05,'Padding',0); plot_isi                (spikes,clusterID);
        
        % Wait for key press
        
        m = 0;
        while ~m
            m   = waitforbuttonpress;
            but = double(get(gcf,'CurrentCharacter'));
        end
        
        % Process key press
        
        if     isequal(but,49) % 1
            spikes.clusters.vars(iClust).manual = 'artifact';
        elseif isequal(but,50) % 2
            spikes.clusters.vars(iClust).manual = 'dubious';
        elseif isequal(but,51) % 3
            spikes.clusters.vars(iClust).manual = 'multi';
        elseif isequal(but,52) % 4
            spikes.clusters.vars(iClust).manual = 'overlap';
        elseif isequal(but,53) % 5
            spikes.clusters.vars(iClust).manual = 'noisy';
        elseif isequal(but,54) % 6
            spikes.clusters.vars(iClust).manual = 'single';
        elseif isequal(but,29) % right arrow
            spikes.clusters.vars(iClust).manual = 'unknown'; % skip 
        elseif isequal(but,28) % left arrow
            iClust = iClust - 2; % go back to previous clusters
        else % stay with current cluster
            iClust = iClust - 1;
        end
    end
    
    iClust = iClust + 1;
    if iClust < 1; iClust = 1; end
end

close(fig);

%------------- END OF CODE --------------