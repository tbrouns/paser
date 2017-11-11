function data = psr_batch_analysis(filenames,loadPath,savePathData,savePathFigs,PLOTTING)

data = []; % Save all data for group analysis

if (nargin < 1 || isempty(filenames))
    filenames = dir('Spikes_*.mat');
    if (size(filenames,1) == 0); return; end
    filenames = char(filenames.name);
end
if (nargin < 2); loadPath     = [];    end
if (nargin < 3); savePathData = [];    end
if (nargin < 4); savePathFigs = [];    end
if (nargin < 5); PLOTTING     = false; end

% Stimulus condition trials
% Stimulus onset trials

close all

fpath = [savePathData '\vars.mat'];
if (exist(fpath,'file') == 2) %% Load data
    load(fpath);
else %% Extract data
    
    params = psr_analysis_parameters();
    
    [filenames,stims] = psr_load_files_session(filenames,loadPath);
    nProbes   = size(filenames,1);
    stimArray = zeros(nProbes,length(stims));
    stims     = unique(stims);
    nStims    = length(stims);
    
    % Arrays to store data in
    
    SpikeCountSingleAll = cell(nStims,nProbes);
    SpikeTimesSingleAll = cell(nStims,nProbes);
    ClusterIDsSingleAll = cell(nStims,nProbes);
    
    SpikeCountPopulationAll = cell(nStims,nProbes);
    SpikeTimesPopulationAll = cell(nStims,nProbes);
    
    for iProbe = 1:nProbes
        
        %% Load file
        
        filename = filenames{iProbe};
        if (isempty(filename)); continue; end
        load([loadPath filename]);
        
        spikes.info.stimulus = metadata.stimulus;
        stimArray(iProbe,:)  = metadata.stimulus;
        
        spikes.info.stimtimes  = metadata.stimtimes;
        spikes.info.trialonset = metadata.trialonset;
        spikes.info.trialonset(end + 1) = spikes.info.dur;
        
        %% Filter clusters
        
        spikes = psr_sst_cluster_filter(spikes,parameters);
        
        %% Extract data for analysis
        
        if (length(metadata.stimulus) > 1)
            for iStim = 1:nStims
                
                %% Take trial data
                
                which = find(spikes.trials == iStim);
                spikesTrl = psr_sst_spike_removal(spikes,which,'keep');
                spikesTrl = psr_sst_spike_removal(spikesTrl,find(spikesTrl.removed),'delete');
                
                spikesTrl.info.stimtimes  = spikesTrl.info.stimtimes (iStim,:);
                spikesTrl.info.dur        = spikesTrl.info.trialonset(iStim + 1) - spikesTrl.info.trialonset(iStim);
                spikesTrl.info.trialonset = spikesTrl.info.trialonset(iStim);
                spikesTrl.info.stimulus   = spikesTrl.info.stimulus  (iStim);
                
                %% Population
                
                spikesPop = spikesTrl;
                type    = [spikesPop.clusters.vars.type];
                tf      = type >= 3; % Take multi-unit and better
                spikesTrl.clusters.vars(~tf) = [];
                ID      = [spikesPop.clusters.vars.id];
                id_all  = [];
                nclusts = length(ID);
                
                for iClust = 1:nclusts
                    id = find(spikesPop.assigns == ID(iClust));
                    id_all = [id_all,id]; %#ok
                end
                
                spikesPop = psr_sst_spike_removal(spikesPop,id_all,'keep');
                spikesPop.assigns(:) = 1; % assigns all spikes to same cluster
                
                spikesPop.clusters.vars    = [];
                spikesPop.clusters.vars.id = 1;
                
                %% Single units
                
                type = [spikesTrl.clusters.vars.type];
                tf   = type >= 5; % Only full single-units
                spikesTrl.clusters.vars(~tf) = [];
                
                %% Extract data
                
                [SpikeArrayPopulation, SpikeTimesPopulation] = psr_get_stimulus_window(spikesPop,params);
                [SpikeArraySingle,     SpikeTimesSingle]     = psr_get_stimulus_window(spikesTrl,params);
                
                ClusterIDsSingle = [spikesTrl.clusters.vars.id];
                
                %                 T_single     = psr_get_spike_times(spikesTrl);
                %                 T_population = psr_get_spike_times(spikesPop);
                
                %% Store data in arrays
                
                SpikeCountSingleAll{iStim,iProbe} = SpikeArraySingle;
                SpikeTimesSingleAll{iStim,iProbe} = SpikeTimesSingle;
                ClusterIDsSingleAll{iStim,iProbe} = ClusterIDsSingle;
                
                SpikeCountPopulationAll{iStim,iProbe} = SpikeArrayPopulation;
                SpikeTimesPopulationAll{iStim,iProbe} = SpikeTimesPopulation;
                
            end
        else % ACTIVE CONDITION
            
        end
        
    end
    
    stimArray = unique(stimArray,'rows');
    if (size(stimArray,1) > 1); disp('Stimulus conditions between probes do not match'); return; end
    params.stims = stimArray;
    
    save([savePathData '\vars.mat'],...
        'params',...
        'SpikeCountSingleAll',...
        'SpikeTimesSingleAll',...
        'ClusterIDsSingleAll',...
        'SpikeCountPopulationAll',...
        'SpikeTimesPopulationAll','-v7.3');
    
end

%%% Array information
% N_all_population = {Number of stimuli conditions x number of probes}
%        sub-array = {1 x 1}
%        sub-array = [number of samples per bin x number of bins x number of trials]
%
%     N_all_single = {Number of stimuli conditions x number of probes}
%        sub-array = {1 x number of clusters}
%        sub-array = [number of samples per bin x number of bins x number of trials]

%% Create folders and save paths

savePathISIPop = [savePathFigs '\ISI\population\'];
savePathISISgl = [savePathFigs '\ISI\single\'];

savePathACFPop = [savePathFigs '\ACF\population\'];
savePathACFSgl = [savePathFigs '\ACF\single\'];

savePathAmpPop = [savePathFigs '\amplitude\population\'];
savePathAmpSgl = [savePathFigs '\amplitude\single\'];

savePathDiffPop = [savePathFigs '\diff\population\'];
savePathDiffSgl = [savePathFigs '\diff\single\'];

savePathPSTHPop = [savePathFigs '\PSTH\population\'];
savePathPSTHSgl = [savePathFigs '\PSTH\single\'];

savePathXCorrSgl = [savePathFigs '\XCorr\single\'];

[~,~,~] = mkdir(savePathISIPop);
[~,~,~] = mkdir(savePathISISgl);

[~,~,~] = mkdir(savePathACFPop);
[~,~,~] = mkdir(savePathACFSgl);

[~,~,~] = mkdir(savePathAmpPop);
[~,~,~] = mkdir(savePathAmpSgl);

[~,~,~] = mkdir(savePathDiffPop);
[~,~,~] = mkdir(savePathDiffSgl);

[~,~,~] = mkdir(savePathPSTHPop);
[~,~,~] = mkdir(savePathPSTHSgl);

[~,~,~] = mkdir(savePathXCorrSgl);

% % ISI
% 
% ISI = ISI_analysis(SpikeTimesPopulationAll,params,savePathISIPop);
% ISI = ISI_analysis(SpikeTimesSingleAll,    params,savePathISISgl,ClusterIDsSingleAll);
% 
% % ACF
% 
% ACF = ACF_analysis(SpikeTimesPopulationAll,params,savePathACFPop);
% ACF = ACF_analysis(SpikeTimesSingleAll,    params,savePathACFSgl,ClusterIDsSingleAll);
% 
% %% Amplitude dependence
% 
% fDiffPopulation = amplitudeAnalysis(SpikeCountPopulationAll, params);
% fDiffSingle     = amplitudeAnalysis(SpikeCountSingleAll,     params);
% 
% amplitudePlotting(fDiffPopulation, params, savePathAmpPop);
% amplitudePlotting(fDiffSingle,     params, savePathAmpSgl, ClusterIDsSingleAll);
% 
% %% Pre- and post-stimulus spiking difference
% 
% SpikeCountDiffPopulation = Diff_analysis(SpikeCountPopulationAll, params);
% SpikeCountDiffSingle     = Diff_analysis(SpikeCountSingleAll,     params);
% 
% Diff_plotting(SpikeCountDiffPopulation, params, savePathDiffPop);
% Diff_plotting(SpikeCountDiffSingle,     params, savePathDiffSgl);

%% PSTH

[SpikeArrayPopulation,SpikeCountPopulation] = PSTH_analysis(SpikeCountPopulationAll, params);
[SpikeArraySingle,    SpikeCountSingle]     = PSTH_analysis(SpikeCountSingleAll,     params);

PSTH_plotting(SpikeCountPopulation, SpikeArrayPopulation, params, savePathPSTHPop);
PSTH_plotting(SpikeCountSingle,     SpikeArraySingle,     params, savePathPSTHSgl, ClusterIDsSingleAll);

%% XCorr

% XCorr = XCorr_analysis(SpikeTimesSingleAll,params,savePathXCorrSgl,ClusterIDsSingleAll);

% JPSTH

% JPSTH_Analysis

%%% TEMP
ACTIVE = true;

if (~ACTIVE)
    
    %% JSPTH calculation
    
    % Merge and average over the stimulus conditions
    
    nSpikesTimeNew = cell(nclusters,nStims);
    
    for iclust = 1:nclusters
        for istim = 2:nStims
            if (size(nSpikesTimeAll{istim},1) >= iclust)
                nSpikesTimeNew{iclust,istim} = nSpikesTimeAll{istim}{iclust};
            end
        end
    end
    
    nstims_total = zeros(nclusters,nclusters);
    JSPTH = cell(nclusters,nclusters);
    
    for istim = 2:nStims % for each stimulus amplitude
        
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
end

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

function SpikeCount = Diff_analysis(SpikeBinCountAll,params)

nStim   = size(SpikeBinCountAll,1);
nProbes = size(SpikeBinCountAll,2);

tMin = params.t_win(1);
tWin = params.t_win(end) - tMin;
nBin = tWin / params.t_bin;
nWin = length(params.diff_win);

% Sum across bin samples in each window

SpikeCountDiffWinAll = [];

for iStim = nStim:-1:1 % iterate backwards for immediate allocation
    for iProbe = nProbes:-1:1
        
        SpikeBinCount = SpikeBinCountAll{iStim,iProbe};
        nClusts = length(SpikeBinCount);
        
        if (nClusts > 0)
            
            for iWin = 1:nWin % Split window in different intervals
                
                SpikeBinCountClustDiff = zeros(size(SpikeBinCount{1},3),nClusts);
                
                % Define sub-window
                iZero  = round((nBin - 1) * ((-tMin) / tWin)) + 1;
                n = ceil((nBin - 1) * ((params.diff_win(iWin) - tMin) / tWin) + 1) - iZero;
                iStart = iZero - n;
                iEnd   = iZero + n;
                
                for iClust = 1:nClusts
                    
                    SpikeBinCountClust = SpikeBinCount{iClust};
                    
                    if (~isempty(SpikeBinCountClust))
                        
                        SpikeBinCountClustPre  = SpikeBinCountClust(:,iStart:iZero,:); % Extract  pre-stimulus window
                        SpikeBinCountClustPost = SpikeBinCountClust(:, iZero:iEnd, :); % Extract post-stimulus window
                        
                        SpikeBinCountClustPre  = sum(SpikeBinCountClustPre,1);
                        SpikeBinCountClustPre  = sum(SpikeBinCountClustPre,2);
                        
                        SpikeBinCountClustPost = sum(SpikeBinCountClustPost,1);
                        SpikeBinCountClustPost = sum(SpikeBinCountClustPost,2);
                        
                        D = SpikeBinCountClustPost - SpikeBinCountClustPre;
                        
                        SpikeBinCountClustDiff(:,iClust) = D(:);
                    end
                end
                
                SpikeCountDiffWinAll{iStim,iProbe,iWin} = mean(SpikeBinCountClustDiff);
            end
        end
    end
end

SpikeCount = SpikeCountDiffWinAll;

end

function Diff_plotting(SpikeDiffAll,params,savePath,clustIDs)

if (nargin < 5); clustIDs = []; end

nProbes = size(SpikeDiffAll,2);
nWin    = size(SpikeDiffAll,3);

stims = params.stims;
[stims,I] = sort(stims);
stimStart = find(stims > 0,1,'first');
stims = stims(stimStart:end);

figure;
colors  = varycolor(nWin);

for iProbe = 1:nProbes
    
    nclusters = size(SpikeDiffAll{end,iProbe,end},2);
    
    for iClust = 1:nclusters
        
        if (~isempty(clustIDs)); clustID = clustIDs{1,iProbe}(iClust);
        else,                    clustID = iClust;
        end
        
        clf;
        hold on;
        
        for iWin = 1:nWin
            
            SpikeDiffTemp = SpikeDiffAll(:,iProbe,iWin);
            SpikeDiffTemp = SpikeDiffTemp(I);
            SpikeDiffTemp = SpikeDiffTemp(stimStart:end);
            nStim         = length(SpikeDiffTemp);
            SpikeDiffStim = zeros(nStim,1);
            for iStim = 1:nStim
                SpikeDiffStim(iStim) = SpikeDiffTemp{iStim}(iClust);
            end
            
            % Plot
            str = ['0-' num2str(params.diff_win(iWin)) ' ms'];
            plot(stims,SpikeDiffStim,'DisplayName',str,...
                'Marker','d',...
                'Color',           colors(iWin,:),...
                'MarkerEdgeColor', colors(iWin,:),...
                'MarkerFaceColor', colors(iWin,:));
            
        end
        
        h = line([stims(1) stims(end)],[0 0],'Color','k','LineStyle',':');
        hasbehavior(h,'legend',false);
        
        yMax = ylim;
        yMax = ceil(max(abs(yMax)));
        ylim([-yMax yMax]);
        
        xlabel('Stimulus amplitude [V]')
        ylabel('Spike count difference (post - pre)')
        
        title(['Probe: ' num2str(iProbe) ', Cluster: ' num2str(clustID)]);
        legend('show','Location','northoutside','Orientation','horizontal');
        export_fig([savePath 'Diff_' num2str(iProbe) '_' num2str(clustID)]);
        
    end
end

end

function [SpikeCount,FiringRate] = PSTH_analysis(SpikeBinCountAll,params)

params.f_win  = [-400 0];
FiringRateAll = calculateFiringRate(SpikeBinCountAll,params);

nStim   = size(SpikeBinCountAll,1);
nProbes = size(SpikeBinCountAll,2);

tMin = params.t_win(1);
tWin = params.t_win(end) - tMin;
tBin = params.t_bin / 1000;
nBin = tWin / params.t_bin;
nWin = length(params.psth_win) - 1;

% Sum across bin samples in each window

FiringRateWinAll = [];
SpikeArrayWinAll = [];

for iStim = nStim:-1:1 % iterate backwards for immediate allocation
    for iProbe = nProbes:-1:1
        
        SpikeBinCount = SpikeBinCountAll{iStim,iProbe};
        nClusts = length(SpikeBinCount);
        
        for iWin = 1:nWin % Split window in different intervals
            
            % Define sub-window
            iStart = round((nBin - 1) * ((params.psth_win(iWin)     - tMin) / tWin)) + 1;
            iEnd   = round((nBin - 1) * ((params.psth_win(iWin + 1) - tMin) / tWin)) + 1;
            
            FiringRateWin = cell(nClusts,1);
            SpikeArrayWin = cell(nClusts,1);
            
            for iClust = 1:nClusts
                
                SpikeBinCountClust = SpikeBinCount{iClust};
                
                if (~isempty(SpikeBinCountClust))
                    
                    SpikeBinCountClust = SpikeBinCountClust(:,iStart:iEnd,:); % Extract sub-window
                    
                    % For raster plot
                    
                    nSamples = size(SpikeBinCountClust,1);
                    nBins    = size(SpikeBinCountClust,2);
                    nTrials  = size(SpikeBinCountClust,3);
                    SpikeArrayWin{iClust} = sparse(reshape(SpikeBinCountClust,nSamples*nBins,nTrials));
                    
                    % For PSTH
                    
                    FiringRateTemp = sum(SpikeBinCountClust,1) / tBin; % Firing rate in each bin
                    FiringRateTemp = permute(FiringRateTemp,[2 3 1]); % Convert to 2D array
                    FiringRateWin{iClust} = FiringRateTemp;
                end
            end
            
            FiringRateWinAll{iStim,iProbe,iWin} = FiringRateWin;
            SpikeArrayWinAll{iStim,iProbe,iWin} = SpikeArrayWin;
        end
    end
end

% Return arrays
SpikeCount = SpikeArrayWinAll;
FiringRate = cell(nStim,nProbes,nWin);

% Calculate difference between pre- and post-stimulus spike number

for iStim = 1:nStim
    for iProbe = 1:nProbes
        for iWin = 1:nWin
            FiringRateWin = FiringRateWinAll{iStim,iProbe,iWin};
            nClusts = length(FiringRateWin);
            if (nClusts > 0)
                FiringRateWinClustAll = zeros(size(FiringRateWin{1},1),nClusts);
                for iClust = 1:nClusts
                    FiringRateWinClust = FiringRateWin{iClust};
                    FiringRateWinClust = mean(FiringRateWinClust,2); % average over all stimulus onsets (trials)
                    FiringRateWinClust = FiringRateWinClust - FiringRateAll{iStim,iProbe}(iClust); % subtract average number of spikes in bin window
                    FiringRateWinClustAll(:,iClust) = FiringRateWinClust;
                end
                FiringRate{iStim,iProbe,iWin} = FiringRateWinClustAll;
            end
        end
    end
end

end

function PSTH_plotting(FiringRateAll,SpikeArrayAll,params,savePath,clustIDs)

if (nargin < 5); clustIDs = []; end

stims  = params.stims;
T      = params.psth_win(1):params.t_bin:params.psth_win(end);
tWin   = T(end) - T(1);
tStart = T(1);
T      = 0.5 * params.t_bin + T(1:end-1);

nStims  = size(FiringRateAll,1);
nProbes = size(FiringRateAll,2);

figure; set(gcf,'position',get(0,'screensize'));

for iProbe = 1:nProbes
    for iStim = 1:nStims
        
        FiringRate = FiringRateAll{iStim,iProbe};
        SpikeArray = SpikeArrayAll{iStim,iProbe};
        
        nClusts = size(FiringRate,2);
        
        if (~isempty(FiringRate))
            
            for iClust = 1:nClusts
                
                clf;
                
                if (~isempty(clustIDs)); clustID = num2str(clustIDs{iStim,iProbe}(iClust));
                else,                    clustID = 'Population';
                end
                
                FiringRateClust  = FiringRate(:,iClust);
                FiringRateSmooth = gaussianSmoothing(FiringRateClust,T,params.sigma);
                
                % PSTH
                
                h(2) = subplot(2,1,2); hold on;
                bar(T,FiringRateClust, 1.0,'FaceColor','b','EdgeColor','none');
                %                 bar(T,SpikeCountClust, 1.0,'FaceColor','b','EdgeColor','none','FaceAlpha',0.25);
                %                 bar(T,SpikeCountSmooth,1.0,'FaceColor','r','EdgeColor','none','FaceAlpha',0.75);
                yMax = ceil(max(abs(FiringRateSmooth)) / 25) * 25;
                if (yMax < 1); yMax = 1; end
                line([0 0],[-yMax yMax],'Color','k','LineStyle','--','LineWidth',1.5);
                xlabel('Time [ms]');
                ylabel('Firing rate [spikes / s]');
                ylim([-yMax yMax]);
                
                % Raster plot
                
                SpikeArrayClust = SpikeArray{iClust};
                
                h(1) = subplot(2,1,1); hold on;
                sWin = size(SpikeArrayClust,1);
                [spiketimes,trial] = find(SpikeArrayClust);
                %                 N = length(trial);
                %                 I = randperm(N);
                %                 if (N > params.Nspikes); I = I(1:params.Nspikes); end
                %                 spiketimes = spiketimes(I);
                %                 trial      = trial(I);
                spiketimes = tWin * (spiketimes / sWin) + tStart;
                scatter(spiketimes,trial,3,'filled','MarkerFaceColor','k');
                
                yMax = ylim; yMax = yMax(2);
                line([0 0],[0 yMax],'Color','k','LineStyle','--','LineWidth',1.0);
                ylabel('Trial #');
                linkaxes([h(1),h(2)],'x')
                psr_plot_stacking(h);
                
                % Save
                
                title(['Probe: ' num2str(iProbe) ', Stimulus: ' num2str(stims(iStim)) 'V, Cluster: ' clustID]);
                export_fig([savePath 'PSTH-P' num2str(iProbe) '-C-' clustID '-V' num2str(stims(iStim))]);
                
            end
        end
    end
end

end

function FiringRateDiff = amplitudeAnalysis(SpikeCountAll,params)

nStim   = size(SpikeCountAll,1);
nProbes = size(SpikeCountAll,2);

tMin = params.t_win(1);
tWin = params.t_win(end) - tMin;
nBin = tWin / params.t_bin;
nWin = length(params.t_array) - 1;

% Mean across all bins and sum across bin samples in each window

FiringRateWin = cell(nStim,nProbes,nWin);

for iStim = nStim:-1:1 % iterate backwards for immediate allocation
    for iProbe = nProbes:-1:1
        
        SpikeCount = SpikeCountAll{iStim,iProbe};
        
        for iWin = 1:nWin % Split main-window in different sub-windows
            
            % Define sub-window
            tStart = params.t_array(iWin);
            tEnd   = params.t_array(iWin + 1);
            iStart = round((nBin - 1) * ((tStart - tMin) / tWin)) + 1;
            iEnd   = round((nBin - 1) * ((tEnd   - tMin) / tWin)) + 1;
            tDur   = params.t_bin * (iEnd - iStart + 1) / 1000; % duration of sub-window [sec]
            
            % Split per cluster
            nClusts = length(SpikeCount);
            FiringRateClusterAll = [];
            for iClust = 1:nClusts
                SpikeCountCluster = SpikeCount{iClust};
                if (~isempty(SpikeCountCluster))
                    SpikeCountCluster = SpikeCountCluster(:,iStart:iEnd,:); % Extract sub-window
                    SpikeCountCluster = sum(SpikeCountCluster,1); % Sum all bin samples for each bin
                    SpikeCountCluster = sum(SpikeCountCluster,2); % Sum over all bins in sub-window
                    FiringRateCluster = SpikeCountCluster / tDur; % Firing rate in sub-window
                    FiringRateClusterAll = [FiringRateClusterAll,FiringRateCluster(:)]; %#ok, Combine results from multiple cluster in one array
                end
            end
            
            FiringRateWin{iStim,iProbe,iWin} = FiringRateClusterAll;
            
        end
    end
end

% Calculate difference between pre- and post-stimulus spike number

FiringRateDiff = cell(size(FiringRateWin));

for iStim = 1:nStim
    for iProbe = 1:nProbes
        for iWin = 2:nWin
            
            FiringRate = FiringRateWin{iStim,iProbe,iWin};
            
            if (~isempty(FiringRate))
                
                % Firing rate difference
                fDiffPre = FiringRate - FiringRateWin{iStim,iProbe,1}; % Relative to pre-stimulus window
                
                % Average over all stimulus onsets (trials)
                FiringRateDiff{iStim,iProbe,iWin} = mean(fDiffPre,1);
                
            end
        end
    end
end

end

function amplitudePlotting(FiringRatesAll,params,savePath,clustIDs)

if (nargin < 4); clustIDs = []; end

nProbes = size(FiringRatesAll,2);
nWin    = size(FiringRatesAll,3);
colors  = varycolor(nWin);

stims     = params.stims;
[stims,I] = sort(stims);
stimStart = find(stims > 0,1,'first');
stims     = stims(stimStart:end);

iWinStart = find(params.t_array >= 0,1,'first');

%% Plotting

figure;
for iProbe = 1:nProbes
    nclusters = size(FiringRatesAll{end,iProbe,end},2);
    for iClust = 1:nclusters
        
        if (~isempty(clustIDs)); clustID = num2str(clustIDs{1,iProbe}(iClust));
        else,                    clustID = 'Population';
        end
        
        clf;
        hold on;
        for iWin = iWinStart:nWin
            
            FiringRateTemp = FiringRatesAll(:,iProbe,iWin);
            FiringRateTemp = FiringRateTemp(I);
            FiringRateTemp = FiringRateTemp(stimStart:end);
            nStim = length(FiringRateTemp);
            FiringRateStims = zeros(nStim,1);
            for iStim = 1:nStim
                FiringRateStims(iStim) = FiringRateTemp{iStim}(iClust);
            end
            
            if (~isempty(FiringRateStims))
                str = [num2str(params.t_array(iWin)) '-' num2str(params.t_array(iWin+1)) ' ms'];
                plot(stims,FiringRateStims,'DisplayName',str,...
                    'Marker','d',...
                    'Color',           colors(iWin,:),...
                    'MarkerEdgeColor', colors(iWin,:),...
                    'MarkerFaceColor', colors(iWin,:));
            end
        end
        
        h = line([stims(1) stims(end)],[0 0],'Color','k','LineStyle',':');
        hasbehavior(h,'legend',false);
        
        xlabel('$\bf{Stimulus \ amplitude \ [V]}$',     'Interpreter','Latex');
        ylabel('$\bf{\Delta Mean \ firing \ rate \ increase \ [spikes / s]}$', 'Interpreter','Latex');
        ymax = ylim;
        ymax = ceil(max(abs(ymax)));
        ylim([-ymax,ymax]);
        title(['Probe: ' num2str(iProbe) ', Cluster: ' clustID]);
        legend('show','Location','northoutside','Orientation','horizontal');
        export_fig([savePath 'AMP_' num2str(iProbe) '_' clustID]);
    end
end

end

function ISI = ISI_analysis(SpikeTimesAll,params,savePath,clustIDs)

if (nargin < 4); clustIDs = []; end

ISI = [];
stims = params.stims;

nStims  = size(SpikeTimesAll,1);
nProbes = size(SpikeTimesAll,2);

fig = figure;
figPos = [0 0 params.fig_wdth params.fig_hght];
set(fig,'Position',figPos);

for iProbe = 1:nProbes
    for iStim = 1:nStims
        
        SpikeTimes = SpikeTimesAll{iStim,iProbe};
        nclusts = size(SpikeTimes,2);
        
        for iClust = 1:nclusts
            
            if (~isempty(clustIDs)); clustID = num2str(clustIDs{iStim,iProbe}(iClust));
            else,                    clustID = 'Population';
            end
            
            [ISI_pre,~]     = ISI_calculation(SpikeTimes{1,iClust},params);
            [ISI_pst,edges] = ISI_calculation(SpikeTimes{2,iClust},params);
            
            ISI_diff = ISI_pst - ISI_pre;
            
            clf;
            t = 0.5 * diff(edges) + edges(1:end-1);
            bar(t,ISI_diff);
            
            yMax = ylim;
            yMax = max(abs(yMax));
            yMax = ceil(100 * yMax) / 100;
            ylim([-yMax yMax]);
            
            xlabel('Time [s]');
            ylabel('Count difference: post - pre');
            title(['Cluster: ' clustID ', Stim: ' num2str(stims(iStim)) 'V']);
            export_fig([savePath 'ISI_Probe_' num2str(iProbe) '_' clustID '_V' num2str(stims(iStim))]);
            
            % N = length(SpikeCounts);
            
        end
    end
end

end

function [ISI,edges] = ISI_calculation(SpikeTimes,params)

ISIs = diff(SpikeTimes);
edges = (0:params.isi_bin:params.isi_max) / 1000;
ISI = histcounts(ISIs,edges, 'Normalization', 'probability');

end

function ACF = ACF_analysis(SpikeTimesAll,params,savePath,clustIDs)

if (nargin < 4); clustIDs = []; end

ACF = [];
stims = params.stims;
nStims  = size(SpikeTimesAll,1);
nProbes = size(SpikeTimesAll,2);

fig = figure;
figPos = [0 0 params.fig_wdth params.fig_hght];
set(fig,'Position',figPos);

for iProbe = 1:nProbes
    for iStim = 1:nStims
        
        SpikeTimes = SpikeTimesAll{iStim,iProbe};
        nclusts = size(SpikeTimes,2);
        
        for iClust = 1:nclusts
            
            if (~isempty(clustIDs)); clustID = num2str(clustIDs{iStim,iProbe}(iClust));
            else,                    clustID = 'Population';
            end
            
            Tmax = SpikeTimes{3,iClust};
            [ACF_pre,~]    = ACF_calculation(SpikeTimes{1,iClust},params,Tmax);
            [ACF_pst,lags] = ACF_calculation(SpikeTimes{2,iClust},params,Tmax);
            
            ACF_diff = ACF_pst - ACF_pre;
            
            clf;
            plot(lags,ACF_diff,'.');
            xlabel('Lag [ms]');
            ylabel('Auto-correlation difference');
            title(['Cluster: ' clustID ', Stim: ' num2str(stims(iStim)) 'V']);
            ylim([-0.1 0.1]);
            export_fig([savePath 'ACF_Probe_' num2str(iProbe) '_' clustID '_V' num2str(stims(iStim))]);
            
            %             N = length(SpikeCounts);
            
        end
    end
end

end

function [ACF,lags] = ACF_calculation(SpikeTimes,params,Tmax)

edges = 0:(params.acf_bin / 1000):Tmax;
SpikeCounts = histcounts(SpikeTimes,edges);
[ACF,lags] = xcorr(SpikeCounts-mean(SpikeCounts),params.acf_max,'coeff');

end

function XCorr = XCorr_analysis(SpikeTimesAll,params,savePath,clustIDs)

XCorr = [];
nStims  = size(SpikeTimesAll,1);
nProbes = size(SpikeTimesAll,2);

stims = params.stims;

fig = figure;
figPos = [0 0 params.fig_wdth params.fig_hght];
set(fig,'Position',figPos);

for iProbe = 1:nProbes
    for iStim = 1:nStims
        
        SpikeTimes = SpikeTimesAll{iStim,iProbe};
        nclusts = size(SpikeTimes,2);
        
        for iClust = 1:nclusts
            for jClust = 1:nclusts
                
                if (iClust == jClust); continue; end
                
                if (~isempty(clustIDs))
                    clustIDx = clustIDs{iStim,iProbe}(iClust);
                    clustIDy = clustIDs{iStim,iProbe}(jClust);
                else
                    clustIDx = iClust;
                    clustIDy = jClust;
                end
                
                Tmax = SpikeTimes{3,iClust};
                [XCorr_pre,~]    = XCorr_calculation(SpikeTimes{1,iClust},SpikeTimes{1,jClust},params,Tmax);
                [XCorr_pst,lags] = XCorr_calculation(SpikeTimes{2,iClust},SpikeTimes{2,jClust},params,Tmax);
                
                XCorr_diff = XCorr_pst - XCorr_pre;
                
                clf;
                plot(lags,XCorr_diff,'.');
                xlabel('Lag [ms]');
                ylabel('Cross-correlation difference');
                title(['Cluster pair: ' num2str(clustIDx) '-' num2str(clustIDy) ', Stim: ' num2str(stims(iStim)) 'V']);
                ylim([-0.1 0.1]);
                export_fig([savePath 'XCorr_' num2str(iProbe) '_C' num2str(clustIDx) '-' num2str(clustIDy) '_' num2str(stims(iStim))]);
                
                %             N = length(SpikeCounts);
            end
        end
    end
end

end

function [XCorr,lags] = XCorr_calculation(SpikeTimes_X,SpikeTimes_Y,params,Tmax)

edges = 0:(params.acf_bin / 1000):Tmax;
SpikeCounts_X = histcounts(SpikeTimes_X,edges);
SpikeCounts_Y = histcounts(SpikeTimes_Y,edges);
SpikeCounts_X = SpikeCounts_X - mean(SpikeCounts_X);
SpikeCounts_Y = SpikeCounts_Y - mean(SpikeCounts_Y);
[XCorr,lags] = xcorr(SpikeCounts_X,SpikeCounts_Y,params.acf_max,'coeff');

end

function x = gaussianSmoothing(x,T,sigma)

n = ceil(sigma / mean(diff(T)));
g = gausswin(n); % create Gaussian smoothing window
g = g / sum(g); % normalize
x = conv(x, g, 'same');

end

function FiringRateAll = calculateFiringRate(SpikeBinCountAll,params)

nStim   = size(SpikeBinCountAll,1);
nProbes = size(SpikeBinCountAll,2);

tMin = params.t_win(1);
tWin = params.t_win(end) - tMin;
nBin = tWin / params.t_bin;

% Sum across bin samples in each window

FiringRateAll = [];

for iStim = nStim:-1:1 % iterate backwards for immediate allocation
    for iProbe = nProbes:-1:1
        
        SpikeBinCount = SpikeBinCountAll{iStim,iProbe};
        nClusts = length(SpikeBinCount);
        
        % Define sub-window
        iStart = round((nBin - 1) * ((params.f_win(1) - tMin) / tWin)) + 1;
        iEnd   = round((nBin - 1) * ((params.f_win(2) - tMin) / tWin)) + 1;
        
        FiringRate = zeros(nClusts,1);
        
        for iClust = 1:nClusts
            
            SpikeBinCountClust = SpikeBinCount{iClust};
            
            if (~isempty(SpikeBinCountClust))
                
                SpikeBinCountClust = SpikeBinCountClust(:,iStart:iEnd,:); % Extract sub-window
                SpikeBinCountClust = SpikeBinCountClust(:);
                                                
                tlength = length(SpikeBinCountClust) / params.Fs;
                nspikes = sum(SpikeBinCountClust);   
                
                FiringRate(iClust) = nspikes / tlength;
            end
        end
        
        FiringRateAll{iStim,iProbe} = FiringRate;
    end
end


end
