max_miss = 0.05; 
max_rpvs = 0.05;

tPre  = 1000;
tPost = 1000;
tbin  =  10; % ms
twin  =  tPre + tPost; 
t     = -tPre + 0.5 * tbin : tbin : tPost - 0.5 * tbin;

stims  = [0 50 60 70 80 90 100 110 120];
nstims = length(stims);

MATfiles  = dir('Spikes_*_Passive_T15_*_1_*.mat');

numfiles  = size(MATfiles,1);
filenames = char(MATfiles.name);
nclusters = 0;

for iFile = 1:numfiles
    
    close all
    
    load(filenames(iFile,:));
        
    clusterIDs = cell2mat({clusters.vars.id});
    clusterMax = max(clusterIDs);
    clustFlags = cell2mat({clusters.vars.flag});
    clustFlags = clustFlags .* strcmp('single',{clusters.vars.unit});
    clustFlags = clustFlags .* (cell2mat({clusters.vars.missing}) < max_miss);
    clustFlags = clustFlags .* (cell2mat({clusters.vars.rpv})     < max_rpvs);
    clusters.vars(~clustFlags) = [];

    %% PSTH
    
    if (iFile == 1);
        fRateStimTimeAll = cell(numfiles,1);
        fRateStimAmpsAll = zeros(length(stims), clusterMax);
    end
    
    I = find(stims == metadata.stimulus);
    [fRateStimTime,fRateStimAmps] = createPSTH(spikes,clusters,clusterMax,tPre,tPost,tbin);
    
    if (~isempty(fRateStimTime) || ~isempty(fRateStimAmps));
        fRateStimTimeAll{iFile} = fRateStimTime;
        fRateStimAmpsAll(I,:)   = fRateStimAmps(1:size(fRateStimAmpsAll,2));
        if size(fRateStimTime,2) > nclusters;
            nclusters = size(fRateStimTime,2);
        end
    end
    
    %% Histogram
    
%     createSpikeHistogram(spikes,clusters);
        
    disp(['Number of files remaining: ' num2str(numfiles - iFile + 1)]);
end

%% Plotting

% Plot spike rate dependence on time from stimulus onset

figure;
hold on
ymax = 2;
for icluster = 1:nclusters
    fRates = [];
    for istim = 2:nstims
        if (size(fRateStimTimeAll{1},2) >= icluster && size(fRateStimTimeAll{istim},2) >= icluster)
            fRates = [fRates,fRateStimTimeAll{istim}(:,icluster) ./ fRateStimTimeAll{1}(:,icluster)];
        end
    end
    mu  = mean(fRates,2);
    sd =  std(fRates,[],2) / sqrt(size(fRates,2));
    if (isempty(mu) || sum(isinf(mu)) > 0 || sum(isnan(mu)) > 0); continue; end
    errorbar(t,mu,sd,'-*')
    if (ceil(max(mu+sd)) > ymax);
        ymax = ceil(max(mu+sd));
    end
end

plot([min(t)-0.5*tbin max(t)+0.5*tbin],[1 1],'--k');
axis([min(t)-0.5*tbin max(t)+0.5*tbin 0 ymax]);
xlabel('Time after stimulus onset [ms]')
ylabel('Relative spike rate increase')

% Plot spike rate dependence on stimulus amplitude

id = ~isnan(fRateStimAmpsAll(1,:));
fRateStimAmpsAll = fRateStimAmpsAll(:,id);
nclusts = size(fRateStimAmpsAll,2);

if (nclusts > 0)

    ymax = 0;

    figure;
    for iclust = 1:nclusts
        fRateStimAmpsAll(:,iclust) = fRateStimAmpsAll(:,iclust) / fRateStimAmpsAll(1,iclust);
        plot(stims(2:end),fRateStimAmpsAll(2:end,iclust),'Marker','*','LineStyle','-');
        hold on
    end

    plot(stims(2:end),nanmean(fRateStimAmpsAll(2:end,:),2),'Marker','*','LineStyle','-','LineWidth',2,'Color','k');

    xlim([min(stims(2:end)) max(stims)]);
    ylim([0 max(max(fRateStimAmpsAll))]);
    xlabel('Stimulus intensity [V]')
    ylabel('Relative spike rate increase')

end
