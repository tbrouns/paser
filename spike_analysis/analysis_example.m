max_miss = 0.05; % Maximum allowable fraction of missing spikes
max_rpvs = 0.05; % Maximum allowable fraction of RPVS

stims  = [0 50 60 70 80 90 100 110 120]; % All stimulus magnitudes
nstims = length(stims);
ntets  = 16;

MATfiles  = dir('Spikes_*.mat'); % Load all MAT files from session in current working directory 
filenames = char(MATfiles.name);

for iTetrode = 1:ntets
    for iStim = 1:nstims
        
        iFile = iStim + (iTetrode - 1) * nstims;
        
        load(filenames(iFile,:));
        
        % Only take complete single unit clusters
        clusterIDs = cell2mat({clusters.vars.id});
        clusterMax = max(clusterIDs);
        clustFlags = cell2mat({clusters.vars.flag});
        clustFlags = clustFlags .* strcmp('single',{clusters.vars.unit});
        clustFlags = clustFlags .* (cell2mat({clusters.vars.missing}) < max_miss);
        clustFlags = clustFlags .* (cell2mat({clusters.vars.rpv})     < max_rpvs);
        clusters.vars(~clustFlags) = [];
        
        I = find(stims == metadata.stimulus); % I: index of "stims" vector
        
        clusterIDs = cell2mat({clusters.vars.id});
        numclusts  = length(clusterIDs);
        
        dataVector = NaN(clusterMax,1); % Vector to contain some data about each cluster
        % Size of 'dataVector' is based on all clusters (single and multi),
        % because we don't know if some clusters are single in every
        % stimulus condition, so we must include every cluster each time.
        % Initialize with NaNs to avoid conflation. 
        
        if (I > 1) % Ignore 0V condition
            % Stimulus onset times in seconds. We have two data vectors here,
            % since we are currently unsure about which of the two intertwined
            % signals indicates the onset times
            stimulusOnset_1 = spikes.artifacts_1;
            stimulusOnset_2 = spikes.artifacts_2;
        end
        
        for iClus = 1:numclusts
            id = (spikes.assigns == clusterIDs(iClus));
            spiketimes = spikes.spiketimes(id); % Only take spike times that are assigned to cluster (given in sec)
            fRate = length(spiketimes) / spikes.info.detect.dur; % Do some basic analysis
            dataVector(iClus) = fRate; % Here we simply store the average firing rate 
        end
    end
end
