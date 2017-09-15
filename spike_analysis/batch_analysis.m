function batch_analysis(MATfiles,loadPath,savePath)

max_miss = 0.05;
max_rpvs = 0.05;

tPre  = 100;
tPost = 300;
tbin  =  10; % ms
t     = -tPre + 0.5 * tbin : tbin : tPost - 0.5 * tbin;
nbins = length(t);

stims  = [0 50 60 70 80 90 100 110 120];
nstims = length(stims);

numfiles  = size(MATfiles,1);
filenames = char(MATfiles.name);
nclusters = 0;

filenamesNew = cell(0,0);
for iFile = 1:numfiles
    filename = filenames(iFile,:);
    load([loadPath filename],'metadata');
    i = metadata.tetrode;
    j = find(stims == metadata.stimulus);
    filenamesNew{i,j} = filename; %#ok
end
filenames = filenamesNew;
ntetrodes = size(filenames,1);

for iTetrode = 1:ntetrodes
    
    close all
    
    for iStim = 1:nstims
        
        filename = filenames{iTetrode,iStim};
        load([loadPath filename]);
        
        clusterIDs = cell2mat({clusters.vars.id});
        clusterMax = max(clusterIDs);
        
        %% Global PSTH 
                        
        if (iStim == 1);
            nSpikesTimeAllGlobal = cell(nstims,1);
        end
        
        [~,~,nSpikesTimeGlobal] = getClusterData(spikes,clusters,clusterMax,tPre,tPost,tbin);
        nSpikesTimeAllGlobal{iStim} = nSpikesTimeGlobal;
        
        %% PSTH
        
        % Only consider single units
        
        clustFlags = cell2mat({clusters.vars.flag});
        clustFlags = clustFlags .* strcmp('single',{clusters.vars.unit});
        clustFlags = clustFlags .* (cell2mat({clusters.vars.missing}) < max_miss);
        clustFlags = clustFlags .* (cell2mat({clusters.vars.rpv})     < max_rpvs);
        clusters.vars(~clustFlags) = [];
                
        if (iStim == 1);
            fRateTimeAll   = cell(nstims,1);
            nSpikesTimeAll = cell(nstims,1);
            fRateAmpsAll   = zeros(length(stims), clusterMax);
        end
        
        I = find(stims == metadata.stimulus);
        [fRateTime,fRateAmps,nSpikesTime] = getClusterData(spikes,clusters,clusterMax,tPre,tPost,tbin);
        
        if (~isempty(fRateTime));
            fRateTimeAll{iStim}   = fRateTime;
            nSpikesTimeAll{iStim} = nSpikesTime;
            fRateAmpsAll(I,1:size(fRateAmps,2)) = fRateAmps; %#ok
            if size(fRateTime,2) > nclusters;
                nclusters = size(fRateTime,2);
            end
        end
        
        %% Histogram
        
        %     createSpikeHistogram(spikes,clusters);
        
    end
    
    %% Plotting
    
    [~,filename,~] = fileparts(filename); % used for saving image
    filename = filename(7:end);
    
    %% Global PSTH 
    figure; set(gcf,'position',get(0,'screensize')); hold on;
    ColorSet = varycolor(nstims-1); % requires Add-on
    for iStim = 2:nstims
        nclusts = size(nSpikesTimeAllGlobal{iStim},1);
        total   = zeros(size(nSpikesTimeAllGlobal{iStim}{1}));
        for iClust = 1:nclusts
            N = nSpikesTimeAllGlobal{iStim}{iClust};
            if (size(N) == size(total));
                total = total + N;
            end
        end
        total = mean(total,2);
        plot(t,total,'LineWidth',1.5,'Color',ColorSet(iStim - 1,:));
        xlabel('Time after stimulus onset [ms]')
        ylabel('Mean number of spikes');
    end
    legend('50','60','70','80','90','100','110','120');
    export_fig([savePath 'GSPTH' filename]);

    %% PSTH
    
    fig1 = figure; set(gcf,'position',get(0,'screensize'));
    hold on
    ymax = 2;
    for icluster = 1:nclusters
        % PSTH
        figure(fig1); hold on;
        fRates = [];
        for istim = 2:nstims
            if (size(fRateTimeAll{1},2) >= icluster && size(fRateTimeAll{istim},2) >= icluster)
                fRates = [fRates,fRateTimeAll{istim}(:,icluster) ./ fRateTimeAll{1}(:,icluster)]; %#ok
            end
        end
        mu = mean(fRates,2);
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
    export_fig([savePath 'SPTH' filename]);

    %% JSPTH calculation
    
    % Merge and average over the stimulus conditions
    
    nSpikesTimeNew = cell(nclusters,nstims);
    
    for iclust = 1:nclusters
        for istim = 2:nstims
            if (size(nSpikesTimeAll{istim},1) >= iclust)
                nSpikesTimeNew{iclust,istim} = nSpikesTimeAll{istim}{iclust};
            end
        end
    end
    
    nartifacts_total = zeros(nclusters,nclusters);
    JSPTH = cell(nclusters,nclusters);
    for istim = 2:nstims
        for iclust = 1:nclusters
            nSpikes_1 = nSpikesTimeNew{iclust,istim};
            if (~isempty(nSpikes_1))
                for jclust = iclust:nclusters % test
                    nSpikes_2 = nSpikesTimeNew{jclust,istim};
                    if (~isempty(nSpikes_2))
                        nartifacts = size(nSpikes_2,2);
                        M = zeros(nbins,nbins,nartifacts);
                        for iartifact = 1:nartifacts % do comparison per stimulus onset
                            N1 = nSpikes_1(:,iartifact);
                            N2 = nSpikes_2(:,iartifact);
                            for ibin = 1:nbins
                                n1 = N1(ibin);
                                for jbin = 1:nbins
                                    n2 = N2(jbin);
                                    M(ibin,jbin,iartifact) = min([n1,n2]);
                                end
                            end
                        end
                        M = sum(M,3);
                        if (~isempty(JSPTH{iclust,jclust})); JSPTH{iclust,jclust} = JSPTH{iclust,jclust} + sum(M,3);
                        else                                 JSPTH{iclust,jclust} = sum(M,3);
                        end
                        nartifacts_total(iclust,jclust) = nartifacts_total(iclust,jclust) + nartifacts;
                    end
                end
            end
        end
    end
    
    for iclust = 1:nclusters
        for jclust = iclust:nclusters
            JSPTH{iclust,jclust} = JSPTH{iclust,jclust} / nartifacts_total(iclust,jclust);
        end
    end
    
    %% Plot JSPTHs
    
    [X,Y] = meshgrid(t,t);
    for iclust = 1:nclusters
        for jclust = iclust:nclusters
            Z = JSPTH{iclust,jclust};
            if (~isempty(Z))
                close all
                figure; set(gcf,'position',get(0,'screensize'));
                surf(X,Y,Z);
                view(2);
                xlabel(['Cluster ' num2str(iclust) ' - Time [s]']);
                ylabel(['Cluster ' num2str(jclust) ' - Time [s]']);
                h = colorbar;
                ylabel(h, 'Mean number of coinciding spikes')
                export_fig([savePath 'JSPTH' filename '_C' num2str(iclust,'%02d') '_' num2str(jclust,'%02d')]);
            end
        end
    end
    
    %% Spike rate dependence on stimulus amplitude
    
    id = ~isnan(fRateAmpsAll(1,:));
    fRateAmpsAll = fRateAmpsAll(:,id);
    nclusts = size(fRateAmpsAll,2);
    
    if (nclusts > 0)
                
        figure; set(gcf,'position',get(0,'screensize')); hold on 
        
        for iclust = 1:nclusts
            fRateAmpsAll(:,iclust) = fRateAmpsAll(:,iclust) / fRateAmpsAll(1,iclust);
            plot(stims(2:end),fRateAmpsAll(2:end,iclust),'Marker','*','LineStyle','-');
            
        end
        
        plot(stims(2:end),nanmean(fRateAmpsAll(2:end,:),2),'Marker','*','LineStyle','-','LineWidth',2,'Color','k');
        
        xlim([min(stims(2:end)) max(stims)]);
        ylim([0 max(max(fRateAmpsAll))]);
        xlabel('Stimulus intensity [V]')
        ylabel('Relative spike rate increase')
        export_fig([savePath 'Stim' filename]);
    end
end

end