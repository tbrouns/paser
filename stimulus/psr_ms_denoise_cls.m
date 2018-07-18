function deletedSpikes = psr_ms_denoise_cls(spikes,metadata,parameters)

VISUALIZE = false;
if (VISUALIZE)
    close all;
    psr_parameters_display;
    parameters.display.metrics = false;
    savePathFig = 'G:\Data\electrophysiology\AVP\data_figures\artifacts\';
    fig = figure; set(gcf,'position',get(0,'screensize'));
end

% Spike removal parameters

Tmax    = sum(spikes.info.dur);
Fs      = spikes.Fs;
Nmax    = floor(Fs * Tmax) + 1;
pBin    = round(Fs * parameters.ms.denoise.cls.tbin / 1000); % Convolution bin
pStm    = round(Fs * parameters.ms.denoise.cls.tstm / 1000); % Stimulus window
pArt    = round(Fs * parameters.ms.denoise.cls.tart / 1000); % Artifact window
pSlp    = round(Fs * parameters.ms.denoise.cls.tslp / 1000); % Slope window
pSpk    = round(Fs * parameters.ms.denoise.cls.tspk / 1000);
nChans  = size(spikes.waveforms,3);
nBlocks = size(metadata.blockonset,1);
blockonsets = [0;cumsum(spikes.info.dur)];
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

deletedSpikes = false(size(spikes.spiketimes)); % spikes to delete

% Convert stimulus timings

stimulusTimes = cell(0);
for iBlock = 1:nBlocks
    if (metadata.stimamps(iBlock) > 0)
        stimulusTimes{end+1} = metadata.stimtimes{iBlock,1} + blockonsets(iBlock);
    end
end
stimulusTimes = cat(1, stimulusTimes{:});
if (isempty(stimulusTimes)); return; end
stimulusTimes = stimulusTimes(:,1);

% Convert to data point index

stimulusPoints = round(Fs * stimulusTimes) + 1;

artTrialIDs = bsxfun(@plus,stimulusPoints,pArt(1):pArt(2));
stmTrialIDs = bsxfun(@plus,stimulusPoints,pStm(1):pStm(2));
% artTrialIDs = artTrialIDs(:);
% stmTrialIDs = stmTrialIDs(:);
artTrialIDs(artTrialIDs < 1)    = 1;
stmTrialIDs(stmTrialIDs < 1)    = 1;
artTrialIDs(artTrialIDs > Nmax) = Nmax;
stmTrialIDs(stmTrialIDs > Nmax) = Nmax;

tBin = parameters.ms.denoise.cls.tbin / 1000;
tArt = parameters.ms.denoise.cls.tart / 1000;
tStm = parameters.ms.denoise.cls.tstm / 1000;
dArt = tArt(2) - tArt(1);
dStm = tStm(2) - tStm(1);
artWindow = [stimulusTimes + tArt(1),stimulusTimes + tArt(2)];
stmWindow = [stimulusTimes + tStm(1),stimulusTimes + tStm(2)];
nTrials = size(stimulusTimes,1);

spikeIDs_art = getSpikeIDs(artWindow,spikes.spiketimes);
spikeIDs_stm = getSpikeIDs(stmWindow,spikes.spiketimes);
spikeIDs_stm = spikeIDs_stm & ~spikeIDs_art;

unitIDs = unique(spikes.assigns_prior);
nUnits  = length(unitIDs);

for iUnit = 1:nUnits
    
    unitID = unitIDs(iUnit);
    
    if (VISUALIZE)
        saveNameFig = [metadata.subject '_' cell2mat(join(metadata.session,'-')) '_P' num2str(metadata.probe,'%02d') '_C' num2str(unitID,'%03d')];
    end
    
    spikeIDs_unit = spikes.assigns_prior == unitID;
    
    spikeIDs_unit_art = spikeIDs_art & spikeIDs_unit;
    spikeIDs_unit_stm = spikeIDs_stm & spikeIDs_unit;
    
    % Ignore clusters that have excessive amount of spikes in stimulus window
    % (likely artifact cluster)
    
    Fr_ex = sum(spikeIDs_unit_stm) / (nTrials * (dStm - dArt)); % Firing rate outside artifact window
    Fr_in = sum(spikeIDs_unit_art) / (nTrials *         dArt);  % Firing rate  inside artifact window
    Fr_ratio = Fr_in / Fr_ex;
    
    Nspikes = length(spikeIDs_unit_art);
    if (Nspikes < parameters.ms.denoise.cls.spkmin); continue; end
    
    if (Fr_ratio >= parameters.ms.denoise.cls.fc_upper)
        
        % Ignore units that have excessive amount of spikes in the artifact window
        % (likely artifact cluster)
        
        deletedSpikes(spikeIDs_unit) = true;
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
        if (VISUALIZE)
            plotWaveforms(fig,spikes,unitID,parameters,spikeIDs_unit_stm,spikeIDs_unit_art,[]);
            suplabel(['F_{ratio} = ' num2str(Fr_ratio)],'t');
            saveStr = [savePathFig saveNameFig '_REM'];
            export_fig(saveStr);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    elseif (Fr_ratio >= parameters.ms.denoise.cls.fc_lower)
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
        if (VISUALIZE)
            plotWaveforms(fig,spikes,unitID,parameters,spikeIDs_unit_stm,spikeIDs_unit_art,[]);
            suplabel(['F_{ratio} = ' num2str(Fr_ratio)],'t');
            saveStr = [savePathFig saveNameFig];
            export_fig(saveStr);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
        % Find excessive amount of spikes in bin
        
        spikePoints = round(Fs * spikes.spiketimes(spikeIDs_unit) + 1);
        spikeSignal = zeros(Nmax,1);
        spikeSignal(spikePoints) = 1;
        spikeSignalArt  = spikeSignal(artTrialIDs);
        spikeSignalArt  = mean(spikeSignalArt,1); % Average spike count over all trials inside artifact window
        spikeSignalConv = conv(spikeSignalArt,ones(1,pBin),'same'); % Number of spikes in bin inside artifact window
        
        nSpikesBin = Fr_ex * tBin; % Number of spikes per bin outside artifact window
        I_art = find(spikeSignalConv > nSpikesBin * parameters.ms.denoise.cls.fr);
        
        if (isempty(I_art)); continue; end
        
        I_diff = find(diff(I_art) > 1);
        range_off = [I_art(I_diff),I_art(end)     ];
        range_on  = [I_art(1),     I_art(I_diff+1)];
        binRanges = [range_on;range_off]';
        binRanges = (binRanges - 1) / Fs;
        
        nBins = size(binRanges,1);
        
        for iBin = 1:nBins
            
            t1 = artWindow(:,1) + binRanges(iBin,1);
            t2 = artWindow(:,1) + binRanges(iBin,2);
            binDur = mean(t2 - t1);
            binWin = [t1,t2];         
                        
            spikeIDs_bin = getSpikeIDs(binWin,spikes.spiketimes(spikeIDs_unit));
            ids = find(spikeIDs_unit);
            spikeIDs_unit_bin = ids(spikeIDs_bin);
            
            Nspikes = length(spikeIDs_unit_bin);
            
            Fr_in = Nspikes / (nTrials * binDur); % Firing rate in bin
            Fr_ratio = Fr_in / Fr_ex;
                        
            if (Nspikes  <  parameters.ms.denoise.cls.spkmin); continue; end
            if (Fr_ratio >= parameters.ms.denoise.cls.fb_upper)
                
                % Calculate mean error
                waveformsEx = spikes.waveforms(spikeIDs_unit_stm,:,:);
                waveformsIn = spikes.waveforms(spikeIDs_unit_bin,:,:);
                errors = getError(waveformsEx,waveformsIn,pSlp,pSpk,parameters);
                errorMean = mean(errors);
                
                if (errorMean > parameters.ms.denoise.cls.dlower)
                    deletedSpikes(spikeIDs_unit_bin) = true;
                end
                                
                %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
                if (VISUALIZE)
                    del = true(size(spikeIDs_unit_bin));
                    plotWaveforms(fig,spikes,unitID,parameters,spikeIDs_unit_stm,spikeIDs_unit_bin,del);
                    suplabel(['F_{ratio} = ' num2str(Fr_ratio) ', Error = ' num2str(errorMean) ', Window = ' num2str(1000 * binDur) ' ms']);
                    saveStr = [savePathFig saveNameFig '_B' num2str(iBin,'%03d')];
                    export_fig(saveStr);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
            elseif (Fr_ratio >= parameters.ms.denoise.cls.fb_lower)
                
                itr = 1;
                SPIKES_DELETED = true;
                spikeIDs = 1:length(spikes.spiketimes);
                spikesTemp = spikes;
                
                spikeIDs_temp_stm = spikeIDs_unit_stm;
                
                while (SPIKES_DELETED)
                    
                    SPIKES_DELETED = false;
                    
                    spikeIDs_temp_unit = spikesTemp.assigns_prior == unitID;
                    spikeIDs_temp_bin  = getSpikeIDs(binWin,spikesTemp.spiketimes(spikeIDs_temp_unit));
                    ids = find(spikeIDs_temp_unit);
                    spikeIDs_unit_bin = ids(spikeIDs_temp_bin);
                                
                    if (length(spikeIDs_unit_bin) < parameters.ms.denoise.cls.spkmin); break; end
                    
                    waveformsEx = spikesTemp.waveforms(spikeIDs_temp_stm,:,:);
                    waveformsIn = spikesTemp.waveforms(spikeIDs_unit_bin,:,:);
                    errors = getError(waveformsEx,waveformsIn,pSlp,pSpk,parameters);
                                        
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
                        IdxDel = spikeIDs_unit_bin(del);
                        
                        if (sum(del) >= parameters.ms.denoise.cls.spkmin) 
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%% Visualize %%%%%%%%%%%%%%%%%%%%%%%%%%
                            if (VISUALIZE)
                                plotWaveforms(fig,spikesTemp,unitID,parameters,spikeIDs_temp_stm,spikeIDs_unit_bin,del);
                                suplabel(['F_{ratio} = ' num2str(Fr_ratio) ', Window = ' num2str(1000 * binDur) ' ms']);
                                saveStr = [savePathFig saveNameFig '_B' num2str(iBin,'%03d') '_I' num2str(itr,'%03d')];
                                export_fig(saveStr);
                                itr = itr + 1;
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % Remove spikes
                            
                            SPIKES_DELETED = true;
                            deletedSpikes(spikeIDs(IdxDel)) = true;
                            spikesTemp = psr_sst_remove_spikes(spikesTemp,IdxDel,'delete');
                            
                            spikeIDs_temp_unit = spikesTemp.assigns_prior == unitID;
                            spikeIDs_temp_bin  = getSpikeIDs(binWin,spikesTemp.spiketimes(spikeIDs_temp_unit));
                            ids = find(spikeIDs_temp_unit);
                            spikeIDs_unit_bin = ids(spikeIDs_temp_bin);

                            errors           (del)     = [];
                            waveformsIn      (del,:,:) = [];
                            spikeIDs         (IdxDel)  = [];
                            spikeIDs_temp_stm(IdxDel)  = [];
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

function plotWaveforms(fig,spikes,unitID,parameters,I_ex,I_in,del)

waveforms = spikes.waveforms;
assigns   = spikes.assigns_prior;

if (islogical(I_ex)); I_ex = find(I_ex); end
if (islogical(I_in)); I_in = find(I_in); end

set(0,'CurrentFigure',fig); clf; hold on;
ax1 = subplot(2,2,1); spikes.waveforms = waveforms(I_ex,      :,:); spikes.assigns = assigns(I_ex);       psr_sst_plot_waveforms(spikes,unitID,parameters);        h = gca; title([h.Title.String ' - All']);
ax2 = subplot(2,2,2); spikes.waveforms = waveforms(I_in,      :,:); spikes.assigns = assigns(I_in);       psr_sst_plot_waveforms(spikes,unitID,parameters);        h = gca; title([h.Title.String ' - Bin']);
ax3 = subplot(2,2,3); spikes.waveforms = waveforms(I_in(~del),:,:); spikes.assigns = assigns(I_in(~del)); psr_sst_plot_waveforms(spikes,unitID,parameters,'line'); h = gca; title([h.Title.String ' - Bin']);
ax4 = subplot(2,2,4); spikes.waveforms = waveforms(I_in( del),:,:); spikes.assigns = assigns(I_in( del)); psr_sst_plot_waveforms(spikes,unitID,parameters,'line'); h = gca; title([h.Title.String ' - Bin']);

axs = [ax1 ax2 ax3 ax4];
linkaxes(axs);

end

function spikeIDs = getSpikeIDs(window,spiketimes)

signIni  = sign(window(:,1) - spiketimes);
signEnd  = sign(window(:,2) - spiketimes);
signMat  = signIni .* signEnd;
spikeIDs = any(signMat < 0,1);

end