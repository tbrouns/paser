function [GPSTH_all,PSTH_all] = psr_batch_analysis(filenames,loadPath,savePathData,savePathFigs,PLOTTING)

close all

if (exist('vars.mat','file') == 2) %% Load data
    load('vars.mat');
else %% Extract data
    if (nargin < 1 || isempty(filenames))
        filenames = dir('Spikes_*.mat');
        if (size(filenames,1) == 0); return; end
        filenames = char(filenames.name);
    end
        
    if (nargin < 2); loadPath     = []; end
    if (nargin < 3); savePathData = []; end
    if (nargin < 4); savePathFigs = []; end
    if (nargin < 5); PLOTTING = false;  end
    
    params = psr_analysis_parameters();
    
    [filenames,stims] = psr_load_files_session(filenames,loadPath);
    ntetrodes = size(filenames,1);
    stims     = unique(stims);
    nstims    = length(stims);
    stimArray = zeros(ntetrodes,nstims);
    
    % Figures
    
    if (PLOTTING)
        fig1 = figure; set(gcf,'position',get(0,'screensize'));
        fig2 = figure; set(gcf,'position',get(0,'screensize'));
        fig3 = figure; set(gcf,'position',get(0,'screensize'));
        fig4 = figure; set(gcf,'position',get(0,'screensize'));
    end
    
    % Arrays to store data in
    
    PSTH_all  = [];
    GPSTH_all = [];
    
    N_all_single          = [];
    F_all_single          = [];
    N_all_population      = [];
    F_all_population      = [];
    
    for iTetrode = 1:ntetrodes
        
        N_tet_single     = cell(nstims,1);
        F_tet_single     = cell(nstims,1);
        N_tet_population = cell(nstims,1);
        F_tet_population = cell(nstims,1);
        
        %% Load file
        
        filename = filenames{iTetrode};
        if (isempty(filename)); continue; end
        load([loadPath filename]);
        
        spikes.info.stimulus  = metadata.stimulus;
        stimArray(iTetrode,:) = metadata.stimulus;
        
        spikes.info.trialonset(end + 1) = spikes.info.detect.dur;
        
        %% Filter clusters
        
        spikes = psr_sst_clusterfilter(spikes,parameters);
        
        %% Extract data for analysis
        
        ACTIVE = 0;
        
        if (length(metadata.stimulus) > 1)
            for iStim = 1:nstims
                which = find(spikes.trials == iStim);
                spikesStim = psr_sst_spike_removal(spikes,which,'keep');
                spikesStim.info.stimtimes  = spikesStim.info.stimtimes(iStim,:);
                spikesStim.info.detect.dur = spikesStim.info.trialonset(iStim + 1) - spikesStim.info.trialonset(iStim);
                spikesStim.info.trialonset = spikesStim.info.trialonset(iStim);
                spikesStim.info.stimulus   = spikesStim.info.stimulus(iStim);
                                
                [N_single,N_population,F_single,F_population] = extractData(spikesStim,params);
                N_tet_single{iStim,1}     = N_single;
                F_tet_single{iStim,1}     = F_single;
                N_tet_population{iStim,1} = N_population;
                F_tet_population{iStim,1} = F_population;
            end
        else % ACTIVE CONDITION
            
        end
        
        N_all_single          = [N_all_single,     N_tet_single];
        F_all_single          = [F_all_single,     F_tet_single];
        N_all_population      = [N_all_population, N_tet_population];
        F_all_population      = [F_all_population, F_tet_population];
        
    end
    
    save('vars.mat','-v7.3');
    
end

stimArray = unique(stimArray,'rows');
if (size(stimArray,1) > 1); disp('Stimulus conditions between probe do not match'); return; end

%% Individual unit analysis

% PSTH

[N_spk_population,N_bin_population] = PSTH_analysis(N_all_population, F_all_population, params);
[N_spk_single,    N_bin_single]     = PSTH_analysis(N_all_single,     F_all_single,     params);

% Amplitude dependence

N_amp_single     = amplitudeAnalysis(N_all_single,     F_all_single,     params);
N_amp_population = amplitudeAnalysis(N_all_population, F_all_population, params);

% For population: sum or average over all probes for each stimulus condition

N_amp_population = cell2mat(N_amp_population);
N_amp_population =     mean(N_amp_population,2); % mean across probes
N_amp_population = num2cell(N_amp_population);

nstims  = size(N_bin_population,1);
nprobes = size(N_bin_population,2);

for istim = 1:nstims
    Nbin = zeros(size(         N_bin_population{istim,1}));
    Nspk = zeros(size(cell2mat(N_spk_population{istim,1})));
    for iprobe = 1:nprobes
        Nbin = Nbin +          N_bin_population{istim,iprobe};
        Nspk = Nspk + cell2mat(N_spk_population{istim,iprobe});
    end
    N_bin_population{istim,1} = Nbin / nprobes;
    N_spk_population{istim,1} = {Nspk};
end

N_bin_population = N_bin_population(:,1);
N_spk_population = N_spk_population(:,1);

%% Pair-wise analysis

XCorrAnalysis()


JPSTH

%% Plotting 

PSTH_plotting(N_bin_population,N_spk_population,params,stimArray(2:end),savePathFigs,'POP');
PSTH_plotting(N_bin_single,    N_spk_single,    params,stimArray(2:end),savePathFigs,'SGL');

amplitudePlotting(N_amp_population,params,stimArray(2:end),savePathFigs,'POP');
amplitudePlotting(N_amp_single,    params,stimArray(2:end),savePathFigs,'SGL');


% %% PSTH analysis
%
% stim_start = 2;
% N_all = zeros(nstims-stim_start,size(t_off,2));
% for iStim = 1:nstims-stim_start
%     N_all(iStim,:) = PSTH_analysis(N_population{iStim+stim_start},t_off);
% end


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
        export_fig([savePathData 'SPTH' filename]);
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
                    export_fig([savePathData 'JSPTH' filename '_C' num2str(iclust,'%02d') '_' num2str(jclust,'%02d')]);
                    
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
            export_fig([savePathData 'Stim' filename]);
        end
    end
end

end

function [N_single,N_population,F_single,F_population] = extractData(spikes,params)

spikes = psr_sst_spike_removal(spikes,find(spikes.removed),'delete');

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

spikes_population = psr_sst_spike_removal(spikes,id_all,'keep');
spikes_population.assigns(:) = 1; % assigns all spikes to same cluster

spikes_population.clusters.vars    = [];
spikes_population.clusters.vars.id = 1;

%% Single units

type = [spikes.clusters.vars.type];
tf   = type >= 5;
spikes.clusters.vars(~tf) = [];

%% Extract data

[N_population,F_population] = psr_analysis_clusterdata(spikes_population,params);
[N_single,    F_single]     = psr_analysis_clusterdata(spikes,params);


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
    z = gaussianSmoothing(z,T,params.sigma);
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

function [NSpikes,NSpikesBin] = PSTH_analysis(N_all,F_all,params)

nstim   = size(N_all,1);
nprobes = size(N_all,2);

twin = params.t_array(end) - params.t_array(1);
nbin = twin / params.t_bin;
nwin = length(params.t_win) - 1;

% Sum across bin samples in each window

N_new_bin = [];
N_new_spk = [];

for istim = nstim:-1:2 % iterate backwards for immediate allocation
    for iprobe = nprobes:-1:1
                
        N = N_all{istim,iprobe};
        nclusts = length(N);
        
        for iwin = 1:nwin % Split window in different intervals
        
            nstart = round(nbin * ((params.t_win(iwin)     - params.t_array(1)) / twin) + 1);
            nend   = round(nbin * ((params.t_win(iwin + 1) - params.t_array(1)) / twin));
                        
            N_temp_bin = cell(nclusts,1);
            N_temp_spk = cell(nclusts,1);
            for iclust = 1:nclusts
                n = N{iclust};
                n = n(:,nstart:nend,:);
                m = reshape(n,size(n,1)*size(n,2),size(n,3));
                N_temp_spk{iclust} = m;
                n = sum(n,1);
                n = permute(n,[2 3 1]);
                N_temp_bin{iclust} = n;
            end
            
            N_new_bin{istim,iprobe,iwin} = N_temp_bin;
            N_new_spk{istim,iprobe,iwin} = N_temp_spk;
        end
    end
end

% Calculate difference between pre- and post-stimulus spike number

for istim = 2:nstim
    for iprobe = 1:nprobes
        f = F_all{istim,iprobe} * (params.t_bin / 1000);
        for iwin = 1:nwin
            N = N_new_bin{istim,iprobe,iwin};
            nclusts = length(N);
            if (nclusts > 0)
                N_new = zeros(size(N{1},1),nclusts);
                for iclust = 1:nclusts
                    n = N{iclust};
                    n = mean(n,2); % average over all stimulus onsets (trials)
                    n = n - f(iclust); % subtract average number of spikes in bin window
                    N_new(:,iclust) = n;
                end
                N_new_bin{istim,iprobe,iwin} = N_new;
            end
        end
    end
end

NSpikesBin = N_new_bin(2:end,:,:);
NSpikes    = N_new_spk(2:end,:,:);

end

function PSTH_plotting(N_bin,N_spk,params,stims,savePath,label)

height = 0.25;
T = params.t_win(1):params.t_bin:params.t_win(end);
twin = T(end) - T(1);
Tstart = T(1);
T = 0.5 * params.t_bin + T(1:end-1);
nstims  = size(N_bin,1);
nprobes = size(N_bin,2);

figure; set(gcf,'position',get(0,'screensize'));

for iProbe = 1:nprobes
    for iStim = 1:nstims
        clf;
        
        spikesAllBin = N_bin{iStim,iProbe};
        spikesAll    = N_spk{iStim,iProbe};
        
        nclusts = size(spikesAllBin,2);
        
        if (~isempty(spikesAllBin))
            
            for iClus = 1:nclusts

                spikesBin = spikesAllBin(:,iClus);
                spikesBinSmoothed = gaussianSmoothing(spikesBin,T,params.sigma);
                
                % PSTH

                h(2) = subplot(2,1,2); hold on;
                bar(T,spikesBin,        1.0,'FaceColor','b','EdgeColor','none','FaceAlpha',0.25);
                bar(T,spikesBinSmoothed,1.0,'FaceColor','r','EdgeColor','none','FaceAlpha',0.75);
                ymax = ylim; ymax = ymax(2);
                line([0 0],[0 ymax],'Color','k','LineStyle','--','LineWidth',1.5);
                xlabel('Time [ms]');

                % Raster plot

                spikes = spikesAll{iClus};

                h(1) = subplot(2,1,1); 
                hold on;

                swin    = size(spikes,1);
                ntrials = size(spikes,2);
                for iTrial = 1:ntrials
                    t = twin * (find(spikes(:,iTrial)) / swin) + Tstart;
                    t = t';
                    plot([t;t],[ones(size(t))*(iTrial-height).';ones(size(t))*(iTrial+height).'],'k','LineWidth',1.5);
                end
                ymax = ylim; ymax = ymax(2);
                line([0 0],[0 ymax],'Color','k','LineStyle','--','LineWidth',1.0);
                ylabel('Trial #');
                linkaxes([h(1),h(2)],'x')
                psr_plot_stacking(h);

                % Save

                title(['Probe: ' num2str(iProbe) ', Stimulus: ' num2str(stims(iStim)) 'V, Cluster: ' num2str(iClus)]);
                export_fig([savePath 'PSTH_' label '_P' num2str(iProbe) '_C' num2str(iClus) '_V' num2str(stims(iStim))]);

            end
        end
    end
end

end

function NSpikes = amplitudeAnalysis(N_all,F_all,params)

nstim   = size(N_all,1);
nprobes = size(N_all,2);

twin = params.t_array(end) - params.t_array(1);
nbin = twin / params.t_bin;
nwin = length(params.t_array) - 1;

% Mean across all bins and sum across bin samples in each window

N_new = [];

for istim = nstim:-1:2 % iterate backwards for immediate allocation
    for iprobe = nprobes:-1:1
                
        N = N_all{istim,iprobe};
                
        for iwin = 1:nwin % Split window in different intervals
        
            nstart = round(nbin * ((params.t_array(iwin)     - params.t_array(1)) / twin) + 1);
            nend   = round(nbin * ((params.t_array(iwin + 1) - params.t_array(1)) / twin));
            
            nclusts = length(N);
            N_temp = [];
            for iclust = 1:nclusts
                n = N{iclust};
                n = n(:,nstart:nend,:);
                n =  sum(n,1); % sum all bin samples
                n = mean(n,2); % average over all bins
                N_temp = [N_temp,n(:)]; %#ok
            end
            
            N_new{istim,iprobe,iwin} = N_temp; 
            
        end
    end
end

% Calculate difference between pre- and post-stimulus spike number

for istim = 2:nstim
    for iprobe = 1:nprobes
        f = F_all{istim,iprobe} * (params.t_bin / 1000);
        for iwin = 2:nwin
            N = N_new{istim,iprobe,iwin};
%             N = N - N_new{istim,iprobe,1}; % relative to pre-stimulus window
            N = N - f;
            N_new{istim,iprobe,iwin} = mean(N,1); % average over all stimulus onsets (trials)
        end
    end
end

NSpikes = N_new(2:end,:,2:end);

end

function amplitudePlotting(N_all,params,stims,savePath,label)

nstim   = size(N_all,1);
nprobes = size(N_all,2);
nwin    = size(N_all,3);
colors  = varycolor(nwin);

[stims,I] = sort(stims);

%% Plotting

figure;
for iprobe = 1:nprobes    
    nclusters = size(N_all{1,iprobe,1},2);
    for iclust = 1:nclusters
        clf;
        hold on;
        for iwin = 1:nwin
            str = [num2str(params.t_array(iwin+1)) '-' num2str(params.t_array(iwin+2)) ' ms'];
            Nstim = zeros(nstim,1);
            for istim = 1:nstim
                Nstim(istim) = N_all{istim,iprobe,iwin}(iclust);
            end
            if (~isempty(Nstim))
                Nstim = Nstim(I,:);
                plot(stims,Nstim,'DisplayName',str,...
                    'Marker','d',...
                    'Color',           colors(iwin,:),...
                    'MarkerEdgeColor', colors(iwin,:),...
                    'MarkerFaceColor', colors(iwin,:));
            end
        end
        xlabel('$\bf{Stimulus \ amplitude \ [V]}$',     'Interpreter','Latex');
        ylabel('$\bf{\Delta Mean \ spike \ increase}$', 'Interpreter','Latex');
        ystep = 1 / params.y_step_amp;
        ymax  = ystep * ylim;
        ymax  = [floor(min(ymax)),ceil(max(ymax))] / ystep;
        ylim(ymax);
        legend('show');
        export_fig([savePath 'AMP_' label '_' num2str(iprobe)]);
    end
end

end

function x = gaussianSmoothing(x,T,sigma)

n = ceil(sigma / mean(diff(T)));
g = gausswin(n); % create Gaussian smoothing window
g = g / sum(g); % normalize
x = conv(x, g, 'same');

end