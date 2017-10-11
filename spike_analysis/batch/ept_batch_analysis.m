function [GPSTH_all,PSTH_all] = ept_batch_analysis(filenames,loadPath,savePath,PLOTTING)

close all

if (nargin < 1 || isempty(filenames))
    filenames = dir('Spikes_*.mat');
    if (size(filenames,1) == 0); return; end
    filenames = char(filenames.name);
end

if (nargin < 2); loadPath = [];    end
if (nargin < 3); savePath = [];    end
if (nargin < 4); PLOTTING = false; end

params = ept_analysis_parameters();
t      = -params.t_pre : params.t_bin : params.t_post;
t_off  = t(1:end-1) + 0.5 * params.t_bin;
nbins  = length(t_off);

[filenames,stims] = ept_load_files_session(filenames,loadPath);
ntetrodes = size(filenames,1);
stims     = sort(unique(stims));
nstims    = length(stims);

% Figures

if (PLOTTING)
    fig1 = figure; set(gcf,'position',get(0,'screensize'));
    fig2 = figure; set(gcf,'position',get(0,'screensize'));
    fig3 = figure; set(gcf,'position',get(0,'screensize'));
    fig4 = figure; set(gcf,'position',get(0,'screensize'));
end

% Arrays to store data in

PSTH_all = [];
GPSTH_all = [];

for iTetrode = 1:ntetrodes
    
    %% Load file
    
    filename = filenames{iTetrode};
    if (isempty(filename)); continue; end
    load([loadPath filename]);
    
    %% Filter clusters
        
    spikes = ept_sst_clusterfilter(spikes,parameters);
    
    %% Extract data for analysis
    
    ACTIVE = 0;
    
    if (length(metadata.stimulus) > 1)
        for iStim = 1:nstims
            I = find(stims == metadata.stimulus(iStim));
            which = find(spikes.trials == iStim);
            spikesStim = ept_sst_spike_removal(spikes,which,'keep');
            spikesStim.info.stimtimes  = spikesStim.info.stimtimes(iStim,:);
            spikesStim.info.trialonset = spikesStim.info.trialonset(iStim);
            [nSpikesTimeGlobal,fRateTime,nSpikesTime,fRateAmps] = extractData(spikesStim,metadata,params);
            
            if (iStim == 1) % Initialize
                nSpikesTimeAllGlobal =  cell(nstims,1);
                fRateTimeAll         =  cell(nstims,1);
                nSpikesTimeAll       =  cell(nstims,1);
                fRateAmpsAll         = zeros(nstims, size(fRateTime,2));
            end
                        
            nSpikesTimeAllGlobal{I} = nSpikesTimeGlobal;
            
            if (~isempty(fRateTime))
                fRateTimeAll{I}   = fRateTime;
                nSpikesTimeAll{I} = nSpikesTime;
                fRateAmpsAll(I,1:size(fRateAmps,2)) = fRateAmps;
            end
            
        end
    else % ACTIVE CONDITION
        
    end
    
    if (isempty(filename)); continue; end
    [~,filename,~] = fileparts(filename); % used for saving image
    filename = filename(7:end);
    
    %% Plotting
    
    if (~ACTIVE)
        
        %% Global PSTH
        
        if (PLOTTING)
            figure(fig1); 
            clf; 
            hold on; 
            ColorSet = varycolor(nstims);
        end
                
        for iStim = 1:nstims
            nclusts = size(nSpikesTimeAllGlobal{iStim},1);
            if (nclusts == 0); continue; end
            jClust = 1;
            total  = [];
            INI    = true;
            for iClust = 1:nclusts
                N = nSpikesTimeAllGlobal{iStim}{iClust};
                
                if (~isempty(N))
                
                    GPSTH{iStim,jClust} = N; %#ok
                    jClust = jClust + 1;
                    
                    if (INI && ~isempty(N))
                        total = zeros(size(N));
                        INI = false;
                    end

                    if (size(N) == size(total))
                        total = total + N;
                    end
                
                end
            end
            
            total = mean(total,2);
            
            if (PLOTTING)
                if (~isempty(total))
                    if (iStim > 1)
                        plot(t_off,total,'LineWidth',1.5,'Color',ColorSet(iStim,:),'DisplayName',num2str(stims(iStim)));
                    else
                        y = mean(total);
                        line([t_off(1) t_off(end)], [y y], 'Color','k','LineStyle',':','DisplayName','0');
                    end
                end
            end
            %             if (max(total) > ymax)
            %                 ymax = ceil(max(total));
            %             end
            
        end
        
        if (~isempty(GPSTH))
            GPSTH_all = [GPSTH_all,GPSTH]; %#ok
        end
        
        if (PLOTTING)
            xlabel('Time after stimulus onset [ms]');
            ylabel('Mean number of spikes');
            if (nstims > 1)
                legend('show','Location','southeast');
            end
            export_fig([savePath 'GSPTH' filename]);
        end
                
        %% PSTH
        
        if (PLOTTING); figure(fig2); clf; hold on; end
        
        nclusters = size(fRateTimeAll{1},2);
        ColorSet = varycolor(nclusters); % requires Add-on
        ymax = 2;
        
        for icluster = 1:nclusters
            fRates = [];
            for istim = 2:nstims
                if (size(fRateTimeAll{1},2) >= icluster && size(fRateTimeAll{istim},2) >= icluster)
                    fRates = [fRates,fRateTimeAll{istim}(:,icluster) ./ fRateTimeAll{1}(:,icluster)]; %#ok
                end
            end
            mu = mean(fRates,2);
            sd =  std(fRates,[],2) / sqrt(size(fRates,2));
            
            if (isempty(mu) || sum(isinf(mu)) > 0 || sum(isnan(mu)) > 0); continue; end
            
            if (PLOTTING)
                errorbar(t_off,mu,sd,'-*','Color',ColorSet(icluster,:));
            end
            
            if (ceil(max(mu+sd)) > ymax); ymax = ceil(max(mu+sd)); end
            
            %% Save
                       
            PSTH_all = [PSTH_all,mu]; %#ok
            
        end
        
        if (PLOTTING)
            plot([min(t_off)-0.5*params.t_bin max(t_off)+0.5*params.t_bin],[1 1],'--k');
            axis([min(t_off)-0.5*params.t_bin max(t_off)+0.5*params.t_bin 0 ymax]);
            xlabel('Time after stimulus onset [ms]')
            ylabel('Relative spike rate increase')
            export_fig([savePath 'SPTH' filename]);
        end
        
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
        
        nstims_total = zeros(nclusters,nclusters);
        JSPTH = cell(nclusters,nclusters);
        
        for istim = 2:nstims % for each stimulus amplitude
            
            for iclust = 1:nclusters % pair-wise comparions between clusters
                
                nSpikes_1 = nSpikesTimeNew{iclust,istim};
                
                if (~isempty(nSpikes_1))
                    
                    for jclust = iclust+1:nclusters
                        
                        nSpikes_2 = nSpikesTimeNew{jclust,istim};
                        
                        if (~isempty(nSpikes_2))
                            
                            nstim_onsets = size(nSpikes_2,2);
                            
                            M = zeros(nbins,nbins,nstim_onsets);
                            
                            for jstim = 1:nstim_onsets % do comparison per stimulus onset
                                
                                N1 = nSpikes_1(:,jstim); %
                                N2 = nSpikes_2(:,jstim);
                                
                                for ibin = 1:nbins
                                    
                                    n1 = N1(ibin);
                                    
                                    for jbin = 1:nbins
                                        
                                        n2 = N2(jbin);
                                        
                                        M(ibin,jbin,jstim) = min([n1,n2]); % number of co-occurences in bin
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            M = sum(M,3); % sum across all stimuli
                            
                            if (~isempty(JSPTH{iclust,jclust})); JSPTH{iclust,jclust} = JSPTH{iclust,jclust} + M;
                            else,                                JSPTH{iclust,jclust} = M;
                            end
                            
                            nstims_total(iclust,jclust) = nstims_total(iclust,jclust) + nstim_onsets;
                        end
                    end
                end
            end
        end
        
        for iclust = 1:nclusters
            for jclust = iclust+1:nclusters
                JSPTH{iclust,jclust} = 100 * JSPTH{iclust,jclust} / nstims_total(iclust,jclust); % percentage
            end
        end
        
        %% Plot JSPTHs
        
        if (PLOTTING)
            
            figure(fig3);
            
            for iclust = 1:nclusters
                for jclust = iclust+1:nclusters
                    
                    Z = JSPTH{iclust,jclust};
                    
                    if (~isempty(Z))
                        
                        % Add padding
                        
                        Z_min = min(min(Z));
                        Z_max = max(max(Z));
                        
                        Z_pad = Z;
                        Z_pad = [Z_pad;Z_min * zeros(size(Z_pad(1,:)))]; %#ok
                        Z_pad = [Z_pad,Z_min * zeros(size(Z_pad(:,1)))]; %#ok
                        
                        clf;
                        
                        % JSPTH matrix
                        
                        [Tx,Ty] = meshgrid(t,t);
                        
                        subplot(3,3,[1,2,4,5]); 
                        surf(Tx,Ty,Z_pad);
                        view(2);
                        xlabel(['Cluster ' num2str(iclust) ' - Time [ms]']);
                        ylabel(['Cluster ' num2str(jclust) ' - Time [ms]']);
                        h = colorbar;
                        caxis([Z_min Z_max]);
                        ylabel(h, 'Mean number of coinciding spikes');
                        ymax = caxis;
                        ymax = ymax(2);
                        
                        %  Correlation averaged over the stimulus duration
                        
                        [Tx,Ty] = meshgrid(t_off,t_off);
                        
                        I = Tx < 0 & Ty < 0;
                        J = Tx > 0 & Ty > 0;
                        
                        tmax = max(max(Tx));
                                                        
                        Z_pre  = reshape(Z(I), max(sum(I)), max(sum(I,2)));
                        Z_post = reshape(Z(J), max(sum(J)), max(sum(J,2)));
                        
                        subplot(3,3,7); 
                        if (~isempty(Z))   
                            bar_JPSTH(Z,params,'g',ymax); 
                            xlim([-tmax,tmax]); 
                            title('Overall');       
                        end
                        
                        subplot(3,3,8);
                        if (~isempty(Z_pre))
                            bar_JPSTH(Z_pre, params,'b',ymax); 
                            xlim([-tmax,tmax]); 
                            title('Pre-stimulus');  
                        end
                        
                        subplot(3,3,9); 
                        if (~isempty(Z_post)) 
                            bar_JPSTH(Z_post,params,'r',ymax); 
                            xlim([-tmax,tmax]); 
                            title('Post-stimulus');
                        end
                                               
                        % Save
                        export_fig([savePath 'JSPTH' filename '_C' num2str(iclust,'%02d') '_' num2str(jclust,'%02d')]);
                        
                    end
                end
            end
        end
        
        %% Spike rate dependence on stimulus amplitude
        
        if (PLOTTING); figure(fig4); clf; hold on; end
        
        id = ~isnan(fRateAmpsAll(1,:));
        fRateAmpsAll = fRateAmpsAll(:,id);
        nclusts = size(fRateAmpsAll,2);
        
        if (nclusts > 0 && nstims > 1)
            
            for iclust = 1:nclusts
                fRateAmpsAll(:,iclust) = fRateAmpsAll(:,iclust) / fRateAmpsAll(1,iclust);
                if (PLOTTING)
                    plot(stims(2:end),fRateAmpsAll(2:end,iclust),'Marker','*','LineStyle','-');
                end
            end
            
            if (PLOTTING)
                plot(stims(2:end),nanmean(fRateAmpsAll(2:end,:),2),'Marker','*','LineStyle','-','LineWidth',2,'Color','k');
                xlim([min(stims(2:end)) max(stims)]);
                ymax = max(max(fRateAmpsAll));
                if (isempty(ymax) || isnan(ymax)); ymax = Inf; end
                ylim([0 ymax]);
                xlabel('Stimulus intensity [V]')
                ylabel('Relative spike rate increase')
                export_fig([savePath 'Stim' filename]);
            end
        end
    end
end

end

function [nSpikesTimeGlobal,fRateTime,nSpikesTime,fRateAmps] = extractData(spikes,metadata,params)

spikes   = ept_sst_spike_removal(spikes,find(spikes.removed),'delete');
clusters = spikes.clusters;

if (metadata.stimulus == 0); control = 1; % check if control trial
else,                        control = 0;
end

clusterIDs = [clusters.vars.id];
clusterMax = max(clusterIDs);

%% Global PSTH

% All spiking units

clustersMulti = clusters;
type = [clusters.vars.type];
tf   = type >= 3;
clustersMulti.vars(~tf) = [];

[~,~,nSpikesTimeGlobal] = ept_analysis_clusterdata(spikes,clustersMulti,control,clusterMax,params);

%% PSTH

% Only consider single units

clustersSingle = clusters;
type = [clusters.vars.type];
tf   = type >= 5;
clustersSingle.vars(~tf) = [];

[fRateTime,fRateAmps,nSpikesTime] = ept_analysis_clusterdata(spikes,clustersSingle,control,clusterMax,params);

end

function bar_JPSTH(Z,params,col,ymax)

Z(Z == 0) = realmax;
z = spdiags(Z); 
z(z == realmax) = 0;
n = size(z,1);
z = fliplr(sum(z));
z = z ./ [1:n,n-1:-1:1];
tlim = (n - 1) * params.t_bin;

T = linspace(-tlim,tlim,length(z));

if (params.smooth)
    n = ceil(params.sigma / mean(diff(T)));
    g = gausswin(n); % create Gaussian smoothing window
    g = g / sum(g); % normalize
    z = conv(z, g, 'same');
end

p = params.pad;

z = z(p:end-p+1);
T = T(p:end-p+1);

bar(T, z, 1.0,'FaceColor',col,'FaceAlpha',0.5);
ylim([0 ymax]);

set(gca,'XTick',-params.t_post:50:params.t_post);
xlabel('Time [ms]')
ylabel('Num. of spikes');

end