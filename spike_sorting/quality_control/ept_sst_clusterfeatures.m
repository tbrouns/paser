function clusters = ept_sst_clusterfeatures(spikes,freq,parameters)

parameters  = ept_parameter_config(); % TEMP
cluster_min = parameters.cluster.size_min;
threshold   = parameters.cluster.thresh * (spikes.info.detect.thresh / spikes.params.detect.thresh);
artifacts   = freq.artifacts;
artifacts   = sparse(getArtifactVector(spikes,artifacts));

% Calculate features of clusters

clusterIDs = unique(spikes.assigns);
nclusters  = length(clusterIDs);

clusters.bhattacharyya = NaN(nclusters,nclusters);
clusters.mahal         = NaN(nclusters,nclusters);

for icluster = 1:nclusters
    
    clusterID = clusterIDs(icluster);
    nspikes   = sum(spikes.assigns == clusterID);
    
    [~,~,~,rpvs] = ss_rpv_contamination (spikes, clusterID); % fraction of RPVs
    [p,~,~,~,~]  = ss_undetected        (spikes, clusterID); % missing spikes due to high threshold
    c            = ss_censored          (spikes, clusterID); % missing spikes due to shadow period
    xc           = crossCorrelation     (spikes, clusterID, parameters); % cross-correlation of residuals
    [k,M,I]      = checkThreshold       (spikes, clusterID, threshold);
    f            = getFiringRate        (spikes, clusterID);
    od           = getArtifactOverlap   (spikes, clusterID, artifacts); % likelihood of spikes being artifacts based on LFP
    snr          = getSignalToNoise     (spikes, clusterID);
    
    % Calculate signal-to-noise ratio
    
    clusters.vars(icluster).id        = clusterID;
    clusters.vars(icluster).rpv       = rpvs / nspikes;
    clusters.vars(icluster).missing   = p;
    clusters.vars(icluster).censored  = c;
    clusters.vars(icluster).nspikes   = nspikes;
    clusters.vars(icluster).xcorr     = xc;
    clusters.vars(icluster).amplitude = M / abs(threshold(I));
    clusters.vars(icluster).chans     = k;
    clusters.vars(icluster).frate     = f;
    clusters.vars(icluster).artifact  = od;
        
    PCA_1 = pca3(spikes,find(spikes.assigns == clusterID));
    for jcluster = (icluster+1):nclusters
        PCA_2 = pca3(spikes,find(spikes.assigns == clusterIDs(jcluster)));
        if (length(PCA_1) > 3 * cluster_min && length(PCA_2) > 3 * cluster_min) % factor of 3 added because pca3 returns three values
            D = mahal(PCA_1,PCA_2);
            clusters.mahal(icluster,jcluster) = mean(D(:));
            [x1,x2,~] = plot_fld(spikes, clusterID, clusterIDs(jcluster), 0);
            clusters.bhattacharyya(icluster,jcluster) = bhattacharyya(x1',x2');
        end
    end
end

% Mixture of drifting t-distribution model for sorting spikes and measuring unit isolation

metrics  = ept_sst_sorting_MDT(spikes,parameters);
ID1      = [metrics.id];
ID2      = [clusters.vars.id];
nClusts  = length(ID1);

for icluster = 1:nClusts
    I = logical(ID2 == ID1(icluster));
    clusters.vars(I).Lratio = metrics(icluster).Lratio;
    clusters.vars(I).IsoDis = metrics(icluster).IsoDis;
    clusters.vars(I).FP_t   = metrics(icluster).FP_t;
    clusters.vars(I).FN_t   = metrics(icluster).FN_t;
    clusters.vars(I).FP_g   = metrics(icluster).FP_g;
    clusters.vars(I).FN_g   = metrics(icluster).FN_g;
end

end

function accept = crossCorrelation(spikes, show, parameters)

PLOTTING = 0;

if (PLOTTING); close all; end

clus      = get_spike_indices(spikes, show);
nchan     = size(spikes.waveforms,3);
nsamples  = size(spikes.waveforms,2);
threshold = parameters.cluster.thresh_xcorr * (spikes.info.detect.thresh / spikes.params.detect.thresh);

%% Calculate cross-correlation between consecutive spike channels

waves = spikes.waveforms(clus,:,:);
waves_mean = squeeze(mean(waves,1));
waves = reshape(permute(waves,[2 1 3]),size(waves,1) * size(waves,2),size(waves,3));

chanIDs = zeros(1,nchan);
waves_max = max(abs(waves_mean));
threshold = abs(threshold); % consider both positive and negative peaks
for ichan = 1:nchan
    chanIDs(ichan) = waves_max(ichan) >= threshold(ichan);
end

if (PLOTTING)
    fig1 = figure; set(gcf,'position',get(0,'screensize'));
    subplot(nchan-1,nchan-1,7);
    plot(1:nsamples,waves_mean);
    xlabel('Sample #');
    ylabel('$\bf{Voltage \ (\mu V)}$','Interpreter','Latex');
    l = legend('Channel 1','Channel 2','Channel 3','Channel 4');
    l.Position(2) = 0.35;
end

maxlag = round(0.75 * nsamples);

N  = (nchan / 2) * (nchan - 1);
XC = zeros(2 * maxlag + 1, N);
accept = zeros(1,N);

n  = 1;
for ichan = 1:nchan
    x = waves(:,ichan);
    for jchan = ichan+1:nchan
        y = waves(:,jchan);
        [XC(:,n),lags] = xcorr(x,y,maxlag,'coeff');
        if (PLOTTING)
            kchan = jchan + (ichan - 1) * nchan;
            kchan = kchan - ceil(kchan / nchan);
            subplot(nchan-1,nchan-1,kchan);
            plot(lags,XC(:,n));
            xlabel('Lag');
            ylabel('Cross correlation');
            title(['Channel ' num2str(ichan) ' vs ' num2str(jchan)]);
            ylim([-1 1]);
            line([0 0],ylim,'Color','k','LineStyle',':');
        end
        accept(n) = chanIDs(ichan) & chanIDs(jchan);
        n = n + 1;
    end
end

if (PLOTTING)
    
    %% Calculate cross-correlations for residuals
    
    memberwaves = spikes.waveforms(clus,:);
    sd          = std(memberwaves);
    sd          = reshape(sd,nsamples,nchan); % residuals for each channel separately
    % sd          = [zeros(size(sd));sd;zeros(size(sd))]; % add zero-padding
    maxlag      = nsamples - 1;
    
    fig2 = figure; set(gcf,'position',get(0,'screensize'));
    st = suptitle('Unbiased'); st.Position(2) = -0.02;
    subplot(nchan-1,nchan-1,7);
    plot(1:nsamples,sd);
    xlabel('Sample #');
    ylabel('Residuals');
    l = legend('Channel 1','Channel 2','Channel 3','Channel 4');
    l.Position(2) = 0.35;
    
    fig3 = figure; set(gcf,'position',get(0,'screensize'));
    st = suptitle('Biased'); st.Position(2) = -0.02;
    subplot(nchan-1,nchan-1,7);
    plot(1:nsamples,sd);
    xlabel('Sample #');
    ylabel('Residuals');
    l = legend('Channel 1','Channel 2','Channel 3','Channel 4');
    l.Position(2) = 0.35;
    
    N  = (nchan / 2) * (nchan - 1);
    XC = zeros(2 * maxlag + 1, N);
    n  = 1;
    for ichan = 1:nchan
        x = sd(:,ichan);
        for jchan = ichan+1:nchan
            kchan = jchan + (ichan - 1) * nchan;
            kchan = kchan - ceil(kchan / nchan);
            
            y = sd(:,jchan);
            
            [XC(:,n),lags] = xcorr(x,y,maxlag,'unbiased');
            figure(fig2);
            subplot(nchan-1,nchan-1,kchan);
            plot(lags,XC(:,n));
            xlabel('Lag');
            ylabel('Cross correlation');
            title(['Channel ' num2str(ichan) ' vs ' num2str(jchan)]);
            line([0 0],ylim,'Color','k','LineStyle',':');
            
            [XC(:,n),lags] = xcorr(x,y,maxlag,'biased');
            figure(fig3);
            subplot(nchan-1,nchan-1,kchan);
            plot(lags,XC(:,n));
            xlabel('Lag');
            ylabel('Cross correlation');
            title(['Channel ' num2str(ichan) ' vs ' num2str(jchan)]);
            line([0 0],ylim,'Color','k','LineStyle',':');
            
            n = n + 1;
        end
    end
    
    figure(fig1); export_fig(['Waveforms_XCorr_Biased_Cluster_' num2str(show)]);
    figure(fig2); export_fig(['Residuals_XCorr_Biased_Cluster_' num2str(show)]);
    figure(fig3); export_fig(['Residuals_XCorr_Unbias_Cluster_' num2str(show)]);
    
end

% Ignore low amplitude channels
XC = XC(:,logical(accept));
ncombi = size(XC,2);

% Find maxima
accept = 1;
lagMax = parameters.cluster.lag_max;
for icombi = 1:ncombi
    xc = XC(:,icombi);
    [pks,loc] = findpeaks(xc);
    [~,I] = max(pks);
    loc = lags(loc(I));
    if (loc > lagMax || loc < -lagMax); 
        accept = 0; 
        break; 
    end
end

%

% figure;
% subplot(2,1,1); hold on;
% bar(lags,xc,1.0);
% line([-maxlag maxlag],[offset offset],'Color','red','LineStyle',':');
% xlim([-maxlag maxlag]);

% xc = max(xc,[],2);
% xc = xc - min(xc);

% subplot(1,2,2);
% bar(lags,xc,1.0);

% for ichan = 1:nchan
%     k = ichan; % number of components in Gaussian mixture model
%
%     OPTIONS = statset('MaxIter',250,'Display','final','TolFun',1e-10);
%     GMModel = gmdistribution.fit(xc_offset,k);
% %     GMModel = fitgmdist(xc_offset,k,'Options',OPTIONS);
%     x = lags(:);
%
% %     wts = GMModel.ComponentProportion;
% %     mu  = GMModel.mu;
% %     sig = sqrt(squeeze(GMModel.Sigma));
%
% %
% %     Ygmfit = zeros(size(x));
% %     for i = 1:ichan
% %         Ygmfit = Ygmfit + wts(i)*normpdf(x, mu(i), sig(i));
% %     end
%
%     % Plot
%
%     for i = 1:ichan
%         plot(x, GMModel.PComponents(i) * pdf('normal', x, GMModel.mu(i), sqrt(GMModel.Sigma(i))),'r')
%     end
% end

end

function [k,M,I] = checkThreshold(spikes, show, threshold)

signThresh   = sign(mean(threshold));
clus         = get_spike_indices(spikes, show);
memberwaves  = spikes.waveforms(clus,:,:);
memberwaves  = signThresh * squeeze(mean(memberwaves,1));
M            = max(memberwaves); % maximum amplitude per channel
k            = find(M > abs(threshold)); % channels that cross threshold with mean amplitude
[M,I]        = max(M); % maximum mean channel amplitude

end

function M = pca3(spikes,which)

waves = spikes.waveforms(which,:);

coeff = pca(waves);

x = spikes.waveforms(which,:) * spikes.info.pca.v(:,1);
y = spikes.waveforms(which,:) * spikes.info.pca.v(:,2);
z = spikes.waveforms(which,:) * spikes.info.pca.v(:,3);

M  = [x,y,z];

end

function f = getFiringRate(spikes,show)

clus       = get_spike_indices(spikes, show);
spiketimes = spikes.spiketimes(clus);
dur        = spikes.info.detect.dur;
f          = length(spiketimes) / dur; % in Hz

end

function artifactVector = getArtifactVector(spikes,artifacts)
    
    Fs        = spikes.params.Fs;
    dur       = spikes.info.detect.dur;
    nlength   = round(Fs * dur);
    artifacts = round(Fs * artifacts);
    artifacts(artifacts < 1) = 1;
    artifacts(artifacts > nlength) = nlength;
    artifactVector = zeros(nlength,1);
    nartifacts = size(artifacts,1);
    for iArtifact = 1:nartifacts
        artifactVector(artifacts(iArtifact,1):artifacts(iArtifact,2)) = 1;
    end
end

function overlap_diff = getArtifactOverlap(spikes,show,artifacts)

% Finds the degree to which spikes overlap with artifacts in LFP signal

Fs         = spikes.params.Fs;
which      = get_spike_indices(spikes, show);
spiketimes = spikes.spiketimes(which);
spiketimes = round(Fs * spiketimes); % in sample number
nlength    = length(artifacts);
nspikes    = length(spiketimes);

n = sum(artifacts);
N = sum(artifacts(spiketimes)); 
overlap_expected = n / nlength; % expected overlap in case of uniform spike distribution
overlap_exact    = N / nspikes; % actual fraction of spikes that overlap with artifact

overlap_diff = overlap_exact / overlap_expected;

end

function snr = getSignalToNoise(spikes,show)

% See: Reliability of Signals From a Chronically Implanted Silicon-Based
% Electrode Array in Non-Human Primate Primary Motor Cortex (2005)

which    = get_spike_indices(spikes, show);
waves    = spikes.waveforms(which,:);
wavesAvg = mean(waves);
A        = max(wavesAvg) - min(wavesAvg);
noise    = bsxfun(@minus,waves,wavesAvg);
noise    = noise(:);
noise    = 2 * std(noise);
snr      = A / noise;

end