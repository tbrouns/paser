function spikes = psr_ms_denoise_cls(spikes,metadata,parameters)

% Temp
psr_parameters_display;
parameters.display.metrics = false;

% Spike removal parameters

Tmax     = sum(spikes.info.dur);
Fs       = spikes.Fs;
Nmax     = floor(Fs * Tmax) + 1;
pBin     = round(Fs * parameters.ms.denoise.cls.tbin / 1000);
pStm     = round(Fs * parameters.ms.denoise.cls.tstm / 1000);
pArt     = round(Fs * parameters.ms.denoise.cls.tart / 1000);
pSlp     = round(Fs * parameters.ms.denoise.cls.tslp / 1000);
pSpk     = round(Fs * parameters.ms.denoise.cls.tspk / 1000);
nChans   = size(spikes.waveforms,3);
nTrials  = size(metadata.trialonset,1);
trialonsets = [0;cumsum(spikes.info.dur)];
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

deletedSpikes = false(size(spikes.spiketimes)); % spikes to delete

% Convert stimulus timings

stimAmps = cell2mat(metadata.stimulus);
% stimAmpsSorted = sort(stimAmps);
% stimMin = stimAmpsSorted(floor(parameters.ms.denoise.cls.ampmin * length(stimAmpsSorted)));

stimulusPoints = cell(0);
for iTrial = 1:nTrials
    if (stimAmps(iTrial) > 0)
        stimulusPoints{end+1} = metadata.stimtimes{iTrial,1} + trialonsets(iTrial);
    end
end
stimulusPoints = cat(1, stimulusPoints{:});
if (isempty(stimulusPoints)); return; end
stimulusPoints = stimulusPoints(:,1);

% Convert to data point index

stimulusPoints = round(Fs * stimulusPoints) + 1;

artTrialIDs = bsxfun(@plus,stimulusPoints,pArt(1):pArt(2));
artTrialIDs(artTrialIDs < 1)    = 1;
artTrialIDs(artTrialIDs > Nmax) = Nmax;

stmTrialIDs = bsxfun(@plus,stimulusPoints,pStm(1):pStm(2));
stmTrialIDs(stmTrialIDs < 1)    = 1;
stmTrialIDs(stmTrialIDs > Nmax) = Nmax;

clustIDs = unique(spikes.assigns_prior);
nClusts  = length(clustIDs);

spikeIDs = 1:length(spikes.spiketimes);
spikesTemp = spikes;

for iClust = 1:nClusts
    
    clustID = clustIDs(iClust);
    
    spikeIDs_clust = find(spikes.assigns_prior == clustID); % Spikes in cluster
    
    [IdxIn,IdxEx] = findIndices(spikes.spiketimes,stmTrialIDs,artTrialIDs,spikeIDs_clust,Nmax,Fs);
    
    % Ignore clusters that have excessive amount of spikes in stimulus window
    % (likely artifact cluster)
    
    Fr_ex = numel(IdxEx) / (numel(stmTrialIDs) - numel(artTrialIDs));
    Fr_in = numel(IdxIn) /                       numel(artTrialIDs);
    
    Fr_ratio = Fr_in / Fr_ex;
    
    Nspikes = length(IdxIn);
    if (Nspikes < parameters.ms.denoise.cls.spkmin); continue; end
    
    if (Fr_ratio >= parameters.ms.denoise.cls.fc_upper)
        
        deletedSpikes(spikeIDs_clust) = true;
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        plotWaveforms(fig,spikes,clustID,parameters,IdxEx,IdxIn,[]);
        suplabel(['F_{ratio} = ' num2str(Fr_ratio)],'t');

        saveStr = [savePath metadata.subject '_' metadata.session{1} '_P' num2str(metadata.probe,'%02d') '_C' num2str(clustID,'%03d') '_REM'];
        export_fig(saveStr);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    elseif (Fr_ratio >= parameters.ms.denoise.cls.fc_lower)
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        plotWaveforms(fig,spikes,clustID,parameters,IdxEx,IdxIn,[]);
        suplabel(['F_{ratio} = ' num2str(Fr_ratio)],'t');

        saveStr = [savePath metadata.subject '_' metadata.session{1} '_P' num2str(metadata.probe,'%02d') '_C' num2str(clustID,'%03d')];
        export_fig(saveStr);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        spikeIDs_clust = find(spikesTemp.assigns_prior == clustID); % Spikes in cluster
        
        [IdxIn,IdxEx] = findIndices(spikesTemp.spiketimes,stmTrialIDs,artTrialIDs,spikeIDs_clust,Nmax,Fs);
        
        % Ignore clusters that have excessive amount of spikes in stimulus window
        % (likely artifact cluster)
        
        Fr_ex = numel(IdxEx) / (numel(stmTrialIDs) - numel(artTrialIDs));
        %     Fr_in = numel(IdxIn) /                       numel(artTrialIDs);
        %
        %     Fr_ratio = Fr_in / Fr_ex;
        %     spikes.clusters.metrics(iClust).stim = Fr_ratio;
        %
        %     if (Fr_ratio > parameters.ms.denoise.cls.ratio)
        %         spikes.clusters.metrics(iClust).quality = 0; % Ignore cluster
        %         continue;
        %     end
        
        % Find excessive amount of spikes in bin
        
        spikePoints = round(Fs * spikesTemp.spiketimes(spikeIDs_clust)) + 1; % Get spike points
        spikePoints(spikePoints > Nmax) = [];
        spikeSignal = zeros(Nmax,1); % Binary vector of spike points
        spikeSignal(spikePoints) = 1;
        spikeSignalArt  = spikeSignal(artTrialIDs);
        spikeSignalArt  = mean(spikeSignalArt,1); % Average spike count over all trials inside artifact window
        spikeSignalConv = conv(spikeSignalArt,ones(1,pBin),'same'); % Number of spikes in bin inside artifact window
        
        nSpikesBin = Fr_ex * pBin; % Number of spikes per bin outside artifact window
        I_art = find(spikeSignalConv > nSpikesBin * parameters.ms.denoise.cls.fr);
        
        if (isempty(I_art)); continue; end
        
        I_diff = find(diff(I_art) > 1);
        range_off = [I_art(I_diff),I_art(end)     ];
        range_on  = [I_art(1),     I_art(I_diff+1)];
        binRanges = [range_on;range_off]';
        
        nBins = size(binRanges,1);
        
        for iBin = 1:nBins
            
            I = (binRanges(iBin,1):binRanges(iBin,2)) + pArt(1);
            binTrialIDs = bsxfun(@plus,stimulusPoints,I);
            binTrialIDs(binTrialIDs < 1)    = 1;
            binTrialIDs(binTrialIDs > Nmax) = Nmax;
            
            binLength = length(I) / Fs;
            [~,IdxEx] = findIndices(spikesTemp.spiketimes,stmTrialIDs,artTrialIDs,spikeIDs_clust,Nmax,Fs);
            IdxIn = spikeIDs_clust(psr_get_spike_ids(spikeSignal,binTrialIDs)); % Inside bin
            Fr_in = numel(IdxIn) / numel(binTrialIDs);
            Fr_ratio = Fr_in / Fr_ex;
            
            Nspikes = length(IdxIn);
            if (Nspikes < parameters.ms.denoise.cls.spkmin); continue; end
    
            if (Fr_ratio >= parameters.ms.denoise.cls.fb_upper)
                
                % Calculate mean error
                waveformsEx = spikes.waveforms(IdxEx,:,:);
                waveformsIn = spikes.waveforms(IdxIn,:,:);
                errors = getError(waveformsEx,waveformsIn,pSlp,pSpk,parameters);
                errorMean = mean(errors);
                
                if (errorMean > parameters.ms.denoise.cls.dlower)
                    deletedSpikes(spikeIDs(IdxIn)) = true;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                del = true(size(IdxIn));
                plotWaveforms(fig,spikesTemp,clustID,parameters,IdxEx,IdxIn,del);
                
                suplabel(['F_{ratio} = ' num2str(Fr_ratio) ', Error = ' num2str(errorMean) ', Window = ' num2str(1000 * binLength) ' ms']);
                saveStr = [savePath metadata.subject '_' metadata.session{1} '_P' num2str(metadata.probe,'%02d') '_C' num2str(clustID,'%03d') '_B' num2str(iBin,'%03d')];
                export_fig(saveStr);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif (Fr_ratio >= parameters.ms.denoise.cls.fb_lower)
                
                itr = 1;
                SPIKES_DELETED = true;
                
                while (SPIKES_DELETED)
                    
                    SPIKES_DELETED = false;
                    
                    [~,IdxEx] = findIndices(spikesTemp.spiketimes,stmTrialIDs,artTrialIDs,spikeIDs_clust,Nmax,Fs);
                    IdxIn = spikeIDs_clust(psr_get_spike_ids(spikeSignal,binTrialIDs)); % Inside bin
                    
                    if (length(IdxIn) < parameters.ms.denoise.cls.spkmin); break; end
                    
                    waveformsEx = spikes.waveforms(IdxEx,:,:);
                    waveformsIn = spikes.waveforms(IdxIn,:,:);
                    errors = getError(waveformsEx,waveformsIn,pSlp,pSpk,parameters);
                    
                    %                 del = error > parameters.ms.denoise.cls.dupper;
                    %
                    %                 artifacts   = waveformsIn(del,:,:);
                    %                 nartifacts  =   size(artifacts,1);
                    %                 artifactMed = median(artifacts,1);
                    %
                    %                 if (nartifacts > thresh); artifactStd = std(artifacts,1);
                    %                 else,                     artifactStd = waveformStd;
                    %                 end
                    %
                    %                 % Calculate difference between artifact median and spike median
                    %
                    %                 d = abs(artifactMed - waveformMed);
                    %                 thresh = parameters.ms.denoise.cls.spkmin;
                    %                 d = d ./ artifactStd;
                    %
                    %                 del = error > parameters.ms.denoise.cls.dupper;
                    %
                    %                 % Classify waveforms inside artifact window
                    %
                    %                 d = calculateError(d, pSpk);
                    %
                    %                 if (d > parameters.ms.denoise.cls.dclus)
                    %
                    %                     diffArt = abs(bsxfun(@minus,artifactMed,waveformsIn));
                    %                     diffSpk = abs(bsxfun(@minus,waveformMed,waveformsIn));
                    %
                    %                     diffArt = mean(diffArt,2);
                    %                     diffSpk = mean(diffSpk,2);
                    %
                    %                     diffArt = sum(weights .* diffArt,3) / weightSum;
                    %                     diffSpk = sum(weights .* diffSpk,3) / weightSum;
                    %
                    %                     del = diffArt < diffSpk;
                    %
                    %                 end
                    
                    [errorMax,errorMaxId] = max(errors);
                    
                    while (~isempty(errors) && errorMax > parameters.ms.denoise.cls.dupper)
                        
                        artifactMax = waveformsIn(errorMaxId,:,:);
                        waveformMed = median(waveformsEx);
                        
                        diffArt = abs(bsxfun(@minus,artifactMax,waveformsIn));
                        diffSpk = abs(bsxfun(@minus,waveformMed,waveformsIn));
                        
                        diffArt = mean(diffArt,2);
                        diffSpk = mean(diffSpk,2);
                        
                        diffArt = mean(diffArt,3);
                        diffSpk = mean(diffSpk,3);
                        
                        del = diffArt < diffSpk;
                        IdxDel = IdxIn(del);
                        
                        if (length(IdxDel) >= parameters.ms.denoise.cls.spkmin)
                            
                            SPIKES_DELETED = true;
                            
                            deletedSpikes(spikeIDs(IdxDel)) = true;
                            
                            spikePoints = round(Fs * spikesTemp.spiketimes(IdxDel)) + 1; % Get spike points
                            spikeSignal(spikePoints) = 0;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            plotWaveforms(fig,spikesTemp,clustID,parameters,IdxEx,IdxIn,del);
                            
                            suplabel(['F_{ratio} = ' num2str(Fr_ratio) ', Window = ' num2str(1000 * binLength) ' ms']);
                            saveStr = [savePath metadata.subject '_' metadata.session{1} '_P' num2str(metadata.probe,'%02d') '_C' num2str(clustID,'%03d') '_B' num2str(iBin,'%03d') '_I' num2str(itr,'%03d')];
                            export_fig(saveStr);
                            itr = itr + 1;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % Remove spikes
                            
                            spikesTemp = psr_sst_remove_spikes(spikesTemp,IdxDel,'delete');
                            spikeIDs_clust = find(spikesTemp.assigns_prior == clustID); % Spikes in cluster
                            
                            errors(del)          = [];
                            waveformsIn(del,:,:) = [];
                            spikeIDs(IdxDel)     = [];
                            
                            %                         IdxEx = find(ismember(spikeIDs,IdxEx));
                            %                         IdxIn = find(ismember(spikeIDs,IdxIn));
                            
                            [~,IdxEx] = findIndices(spikesTemp.spiketimes,stmTrialIDs,artTrialIDs,spikeIDs_clust,Nmax,Fs);
                            IdxIn = spikeIDs_clust(psr_get_spike_ids(spikeSignal,binTrialIDs)); % Inside bin
                            
                        else
                            
                            errors(errorMaxId) = 0; % Ignore error for next iterations
                            
                        end
                        
                        [errorMax,errorMaxId] = max(errors);
                        
                    end
                end
            end
        end
    end
end

spikes = psr_sst_remove_spikes(spikes,deletedSpikes,'delete');

end

function I = calculateRange(nPoints,n)

nPointsHalf = round(0.5 * nPoints);
I1 = nPointsHalf + n(1);
I2 = nPointsHalf + n(2);
i1 = I1; if (i1 < 1);       i1 = 1;       end
i2 = I2; if (i2 > nPoints); i2 = nPoints; end

I = i1:i2;

end

function [I_in,I_ex] = findIndices(spiketimes,stmTrialIDs,artTrialIDs,spikeIDs,Nmax,Fs)

spikePoints = round(Fs * spiketimes(spikeIDs)) + 1; % Get spike points
spikePoints(spikePoints > Nmax) = [];

spikeIDs_stm = psr_get_spike_ids(spikePoints,stmTrialIDs); % Spikes inside stimulus window
spikeIDs_art = psr_get_spike_ids(spikePoints,artTrialIDs); % Spikes inside artifact window

I_ex = spikeIDs(~spikeIDs_art & spikeIDs_stm); % Outside of artifact window, but still inside stimulus window
I_in = spikeIDs( spikeIDs_art);                %  Inside of artifact window

end

function errors = getError(waveformsEx,waveformsIn,pSlp,pSpk,parameters)

    % Calculate slope

    diffWavesEx = waveformsEx(:,1+pSlp:end,:) - waveformsEx(:,1:end-pSlp,:);
    diffWavesIn = waveformsIn(:,1+pSlp:end,:) - waveformsIn(:,1:end-pSlp,:);

    % Calculate average and standard deviation

    waveformMed = median(waveformsEx);
    diffWaveMed = median(diffWavesEx);

    waveformStd = std(waveformsEx); 
    diffWaveStd = std(diffWavesEx);

    % Normalize by standard deviation

    waveformMedNorm = waveformMed ./ waveformStd;
    diffWaveMedNorm = diffWaveMed ./ diffWaveStd;

    if (size(waveformsIn,1) >= parameters.ms.denoise.cls.spkmin)
        waveformStd = std(waveformsIn); 
        diffWaveStd = std(diffWavesIn);
    end
    
    waveformsNorm = bsxfun(@rdivide,waveformsIn,waveformStd);
    diffWavesNorm = bsxfun(@rdivide,diffWavesIn,diffWaveStd);

    % Calculate error

    errorWave  = bsxfun(@minus,waveformMedNorm,waveformsNorm);
    errorSlope = bsxfun(@minus,diffWaveMedNorm,diffWavesNorm);

    errorWave  = calculateError(errorWave, pSpk);
    errorSlope = calculateError(errorSlope,pSpk);
    errors = mean([errorWave,errorSlope],2);
                    
end

function error = calculateError(x,n)

error  = abs(x);
error  = sort(error,2,'descend');
error  = error(:,1:n,:);
error  = mean(error,   2);
error  =  max(error,[],3);

end

function plotWaveforms(fig,spikes,clustID,parameters,I_ex,I_in,del)

waveforms = spikes.waveforms;
assigns   = spikes.assigns_prior;

set(0,'CurrentFigure',fig); clf; hold on;
ax1 = subplot(2,2,1); spikes.waveforms = waveforms(I_ex,      :,:); spikes.assigns = assigns(I_ex);       psr_sst_plot_waveforms(spikes,clustID,parameters);        h = gca; title([h.Title.String ' - All']);
ax2 = subplot(2,2,2); spikes.waveforms = waveforms(I_in,      :,:); spikes.assigns = assigns(I_in);       psr_sst_plot_waveforms(spikes,clustID,parameters);        h = gca; title([h.Title.String ' - Bin']);
ax3 = subplot(2,2,3); spikes.waveforms = waveforms(I_in(~del),:,:); spikes.assigns = assigns(I_in(~del)); psr_sst_plot_waveforms(spikes,clustID,parameters,'line'); h = gca; title([h.Title.String ' - Bin']);
ax4 = subplot(2,2,4); spikes.waveforms = waveforms(I_in( del),:,:); spikes.assigns = assigns(I_in( del)); psr_sst_plot_waveforms(spikes,clustID,parameters,'line'); h = gca; title([h.Title.String ' - Bin']);

axs = [ax1 ax2 ax3 ax4];
linkaxes(axs);

end