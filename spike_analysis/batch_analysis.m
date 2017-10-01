function batch_analysis(filenames,loadPath,savePath)

close all

if (nargin < 1);
    filenames = dir('*.mat');
    if (size(filenames,1) == 0); return; end
    filenames = char(filenames.name);
end

if (nargin < 2); loadPath = []; end
if (nargin < 3); savePath = []; end

tPre  = 100;
tPost = 300;
tbin  =  10; % ms
t     = -tPre + 0.5 * tbin : tbin : tPost - 0.5 * tbin;
nbins = length(t);

sigma = 50; % for smoothing (ms)

nclusters = 0;

[filenames,stims] = ept_load_files_session(filenames,loadPath);
ntetrodes = size(filenames,1);
nstims    = size(filenames,2);

% Figures

fig1 = figure; set(gcf,'position',get(0,'screensize'));
fig2 = figure; set(gcf,'position',get(0,'screensize'));
fig3 = figure; set(gcf,'position',get(0,'screensize'));
fig4 = figure; set(gcf,'position',get(0,'screensize'));

for iTetrode = 1:ntetrodes
    
    if (nstims > 1); % Passive condition
        
        for iStim = 1:nstims
            
            filename = filenames{iTetrode,iStim};
            if (isempty(filename)); continue; end
            load([loadPath filename]);
                        
            parameters = ept_parameter_config(); % TEMP
            
            spikes = ept_sst_spike_removal(spikes,spikes.removed,'delete');
            
            if (metadata.stimulus == 0); control = 1; % check if control trial
            else                         control = 0;
            end
            
            clusterIDs = cell2mat({clusters.vars.id});
            clusterMax = max(clusterIDs);
            
            %% Global PSTH
            
            if (iStim == 1);  nSpikesTimeAllGlobal = cell(nstims,1);  end
            [~,~,nSpikesTimeGlobal] = getClusterData(spikes,clusters,control,clusterMax,tPre,tPost,tbin);
            nSpikesTimeAllGlobal{iStim} = nSpikesTimeGlobal;
            
            %% PSTH
            
            % Only consider single units
            
            rpv = [clusters.vars.rpv] ./ [clusters.vars.nspikes]; % TEMP
            
            nclusts = size([clusters.vars.id],1);
            id = ones(nclusts,1);
            id = id .* ([clusters.vars.missing] < parameters.cluster.upper_miss);
            id = id .* (rpv                     < parameters.cluster.upper_rpv);
            id = logical(id);
            clusters.vars(~id) = [];
            
            if (iStim == 1);
                fRateTimeAll   =  cell(nstims,1);
                nSpikesTimeAll =  cell(nstims,1);
                fRateAmpsAll   = zeros(nstims, clusterMax);
            end
            
            I = find(stims == metadata.stimulus);
            [fRateTime,fRateAmps,nSpikesTime] = getClusterData(spikes,clusters,control,clusterMax,tPre,tPost,tbin);
            
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
        
    else % Active condition
        
        
    end
    
    %% Plotting
    
    if (isempty(filename)); continue; end
    [~,filename,~] = fileparts(filename); % used for saving image
    filename = filename(7:end);
    
    if (nstims > 1) % Passive
        
        %% Global PSTH
        figure(fig1); clf; hold on;
        ymax = 1;
        ColorSet = varycolor(nstims); % requires Add-on
        for iStim = 1:nstims
            nclusts = size(nSpikesTimeAllGlobal{iStim},1);
            if (nclusts == 0); continue; end
            total   = zeros(size(nSpikesTimeAllGlobal{iStim}{1}));
            for iClust = 1:nclusts
                N = nSpikesTimeAllGlobal{iStim}{iClust};
                if (size(N) == size(total));
                    total = total + N;
                end
            end
            
            total = mean(total,2);
            
            if (iStim > 1)
                plot(t,total,'LineWidth',1.5,'Color',ColorSet(iStim,:));
            else
                y = mean(total);
                line([t(1) t(end)], [y y], 'Color','k','LineStyle',':');
            end     
            
            if (max(total) > ymax)
               ymax = ceil(max(total)); 
            end
            
        end
        
        xlabel('Time after stimulus onset [ms]');
        ylabel('Mean number of spikes');
        ylim([0 ymax]);
        
        if (nstims > 1); 
            legend('0','50','60','70','80','90','100','110','120','Location','southeast'); 
        end
        export_fig([savePath 'GSPTH' filename]);
        
        %% PSTH
        figure(fig2); clf; hold on
        ColorSet = varycolor(nclusters); % requires Add-on
        ymax = 2;
        for icluster = 1:nclusters
            % PSTH
            fRates = [];
            for istim = 2:nstims
                if (size(fRateTimeAll{1},2) >= icluster && size(fRateTimeAll{istim},2) >= icluster)
                    fRates = [fRates,fRateTimeAll{istim}(:,icluster) ./ fRateTimeAll{1}(:,icluster)]; %#ok
                end
            end
            mu = mean(fRates,2);
            sd =  std(fRates,[],2) / sqrt(size(fRates,2));
            
            if (isempty(mu) || sum(isinf(mu)) > 0 || sum(isnan(mu)) > 0); continue; end
            
            errorbar(t,mu,sd,'-*','Color',ColorSet(icluster,:));
            
            if (ceil(max(mu+sd)) > ymax); ymax = ceil(max(mu+sd)); end
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
        figure(fig3);
        [X,Y] = meshgrid(t,t);
        for iclust = 1:nclusters
            for jclust = iclust:nclusters
                Z = JSPTH{iclust,jclust};
                if (~isempty(Z))
                    clf;
                    subplot(1,2,1); % JSPTH matrix
                    surf(X,Y,Z);
                    view(2);
                    xlabel(['Cluster ' num2str(iclust) ' - Time [s]']);
                    ylabel(['Cluster ' num2str(jclust) ' - Time [s]']);
                    h = colorbar;
                    ylabel(h, 'Mean number of coinciding spikes');
                    subplot(1,2,2); %  Correlation averaged over the stimulus duration
                    Z(Z == 0) = realmax; z = spdiags(Z); z(z == realmax) = 0;
                    n = size(z,1);
                    z = fliplr(sum(z));
                    z = z ./ [1:n,n-1:-1:1];
                    T = linspace(-t(end),t(end),length(z));
                    
                    n = ceil(sigma / mean(diff(T)));
                    g = gausswin(n); % create Gaussian smoothing window
                    g = g / sum(g); % normalize
                    zs = conv(z, g, 'same');
%                     zs = smooth(T,z,20);

                    T  =  T(5:end-4);
                    z  =  z(5:end-4);
                    zs = zs(5:end-4);
                    
%                     hold on;
%                     bar(T,z, 1.0,'FaceColor','b','FaceAlpha',0.5); 
                    bar(T,zs,1.0,'FaceColor','r','FaceAlpha',0.5);
                    
                    % Save
                    export_fig([savePath 'JSPTH' filename '_C' num2str(iclust,'%02d') '_' num2str(jclust,'%02d')]);
                end
            end
        end
        
        %% Spike rate dependence on stimulus amplitude
        
        id = ~isnan(fRateAmpsAll(1,:));
        fRateAmpsAll = fRateAmpsAll(:,id);
        nclusts = size(fRateAmpsAll,2);
        
        if (nclusts > 0 && nstims > 1)
            
            figure(fig4); clf; hold on
            
            for iclust = 1:nclusts
                fRateAmpsAll(:,iclust) = fRateAmpsAll(:,iclust) / fRateAmpsAll(1,iclust);
                plot(stims(2:end),fRateAmpsAll(2:end,iclust),'Marker','*','LineStyle','-');
            end
            
            plot(stims(2:end),nanmean(fRateAmpsAll(2:end,:),2),'Marker','*','LineStyle','-','LineWidth',2,'Color','k');
            
            xlim([min(stims(2:end)) max(stims)]);
            ymax = max(max(fRateAmpsAll));
            if (isempty(ymax) || isnan(ymax)); ymax = Inf; end
            ylim([0 ymax]);
            xlabel('Stimulus intensity [V]')
            ylabel('Relative spike rate increase')
            export_fig([savePath 'Stim' filename]);
        end
        
    else % Active condition
               
        
    end
end

end