function [GPSTH_all,PSTH_all] = ept_batch_analysis(filenames,loadPath,savePath,PLOTTING)


% load('vars.mat');
skip = false; % TEMP

close all

if (~skip)
    
    if (nargin < 1 || isempty(filenames))
        filenames = dir('Spikes_*.mat');
        if (size(filenames,1) == 0); return; end
        filenames = char(filenames.name);
    end
    
    if (nargin < 2); loadPath = [];    end
    if (nargin < 3); savePath = [];    end
    if (nargin < 4); PLOTTING = false; end
    
    params = ept_analysis_parameters();
    t      = params.t_pre : params.t_bin : params.t_post;
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
    
    N_single_all          = [];
    N_population_all      = [];
    N_population_stim_all = [];
    
    for iTetrode = 1:ntetrodes
        
        N_single_tet     = [];
        N_population_tet = [];
        N_population_stim_tet = [];
        
        %% Load file
        
        filename = filenames{iTetrode};
        if (isempty(filename)); continue; end
        load([loadPath filename]);
        
        spikes.info.stimulus = metadata.stimulus; % TEMP
        
        %% Filter clusters
        
        spikes = ept_sst_clusterfilter(spikes,parameters);
        
        %% Extract data for analysis
        
        ACTIVE = 0;
        
        if (length(metadata.stimulus) > 1)
            for iStim = 1:nstims
                which = find(spikes.trials == iStim);
                spikesStim = ept_sst_spike_removal(spikes,which,'keep');
                spikesStim.info.stimtimes  = spikesStim.info.stimtimes(iStim,:);
                spikesStim.info.trialonset = spikesStim.info.trialonset(iStim);
                spikesStim.info.stimulus   = spikesStim.info.stimulus(iStim);
                [N_single,N_population,N_population_stim] = extractData(spikesStim,params);
                
                nclusts = length(N_single);
                
                if (nclusts > 0)
                    N_single_tet{iStim,nclusts} = []; % allocate space
                    N_single_tet(iStim,:)       = N_single;
                end
                
                if (~isempty(cell2mat(N_population)))
                    N_population_tet{iStim,1} = [];
                    N_population_tet(iStim,1) = N_population;
                end
                
                for i = 1:size(N_population_stim,1)
                    N_population_stim_tet{iStim,1,i} = N_population_stim{i};
                end
            end
        else % ACTIVE CONDITION
            
        end
        
        N_single_all          = [N_single_all,          N_single_tet];
        N_population_all      = [N_population_all,      N_population_tet];
        N_population_stim_all = [N_population_stim_all, N_population_stim_tet];
        
    end
    
end

%% Combine population data and normalize

% PSTH
N_population = combinePopulation(N_population_all,nstims);

% Amplitude dependence

N = amplitudeAnalysis(N_population_stim_all,stims);

% nwin = size(N_population_stim_all,3);
% N_population_stim = cell(nwin,1);
% for iwin = 1:nwin
%     N_population_stim{iwin} = combinePopulation(N_population_stim_all(:,:,iwin),nstims);
% end

%% PSTH analysis

stim_start = 2;
N_all = zeros(nstims-stim_start,size(t_off,2));
for iStim = 1:nstims-stim_start
    N_all(iStim,:) = PSTH_analysis(N_population{iStim+stim_start},t_off);
end

%% Amplitude analysis

N_amps = zeros(nwin,nstims);
for iwin = 1:nwin
    Nwin = N_population_stim{iwin};
    for istim = 1:nstims
        N = Nwin{istim};
        N =  sum(N,1);
        N =  sum(N,2);
        N = mean(N,3);
        N_amps(iwin,istim) = N;
    end
end

%% Plotting
figure; hold on;

for iwin = 1:nwin
    plot(stims(2:end),N_amps(iwin,2:end));
end

if (~ACTIVE)
    
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

function [N_single,N_pop,N_pop_stim] = extractData(spikes,params)

spikes = ept_sst_spike_removal(spikes,find(spikes.removed),'delete');

%% Population

type    = [spikes.clusters.vars.type];
tf      = type >= 3;
ID      = spikes.clusters.vars(tf).id;
id_all  = [];
nclusts = length(ID);

for iClust = 1:nclusts
    id = find(spikes.assigns == ID(iClust));
    id_all = [id_all;id]; %#ok
end

spikes_population = ept_sst_spike_removal(spikes,id_all,'keep');
spikes_population.assigns(:) = 1; % assigns all spikes to same cluster

spikes_population.clusters.vars    = [];
spikes_population.clusters.vars.id = 1;

T = zeros(3,1);
T(1) = params.t_pre;
T(2) = params.t_post;
T(3) = params.t_bin;

N_pop = ept_analysis_clusterdata(spikes_population,T);

%% Stimulus amplitude dependence

n    = length(params.t_array);
t    = zeros(3,1);
t(3) = params.t_bin;
N_pop_stim = cell(n-1,1);

for i = 1:n-1
    t(1) = params.t_array(i);
    t(2) = params.t_array(i + 1);
    N_pop_stim{i} = cell2mat(ept_analysis_clusterdata(spikes_population,t));
end

%% Single units

type = [spikes.clusters.vars.type];
tf   = type >= 5;
spikes.clusters.vars(~tf) = [];

N_single = ept_analysis_clusterdata(spikes,T);

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

function N_new = combinePopulation(N_all,nstims)

nprobes = size(N_all,2);
N_new = cell(nstims,1);
for iStim = 1:nstims
    for iProbe = 1:nprobes
        N = N_all{iStim,iProbe};
        if (iProbe == 1); N_new{iStim} = zeros(size(N)); end
        N_new{iStim} = N_new{iStim} + N ./ N_all{1,iProbe};
    end
end

end

function N = PSTH_analysis(N,t)

N = sum(N);
N = mean(N,3);
N = N - mean(N);

% figure;
% bar(t,N);
% ylim([floor(min(N)) ceil(max(N))]);

end

function N = amplitudeAnalysis(N_all,stims)

nstim  = size(N_all,1);
nprobe = size(N_all,2);
nwin   = size(N_all,3);

% Mean across all bins and all bin samples in each window

for istim = 1:nstim
    for iprobe = 1:nprobe
        for iwin = 1:nwin
            N = N_all{istim,iprobe,iwin};
            N = mean(N,1);
            N = mean(N,2);
            N_all{istim,iprobe,iwin} = N(:);
        end
    end
end

% Calculate difference between pre- and post-stimulus spike number

N_new = zeros(nstim,nprobe,nwin);

for istim = 1:nstim
    for iprobe = 1:nprobe
        for iwin = 2:nwin
            N = N_all{istim,iprobe,iwin};
            N = N - N_all{istim,iprobe,1}; % relative to pre-stimulus window
            N_new(istim,iprobe,iwin) = mean(N);
        end
    end
end

N_new = mean(N_new,2);

figure; hold on;
for iwin = 2:nwin
   plot(stims(2:end),N_new(2:end,:,iwin));
end

end