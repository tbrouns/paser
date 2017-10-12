function SpikeBinTrials = ept_analysis_clusterdata(spikes,T)

tPre  = T(1);
tPost = T(2);
tBin  = T(3); 

tWin  = tPost - tPre; 
Tmax  = spikes.info.detect.dur; 
Fs    = spikes.params.Fs; 
Nmax  = floor(Fs * Tmax);
sPre  = Fs * (tPre  / 1000); % pre- stimulus window
sPost = Fs * (tPost / 1000); % post-stimulus window
sBin  = Fs * (tBin  / 1000); % bin size in stimulus window
nBins = tWin / tBin;

clusterIDs = [spikes.clusters.vars.id];
numclusts  = length(clusterIDs);
SpikeBinTrials = cell(1,numclusts);

% Extract stimulus windows

stimulusAmp   = spikes.info.stimulus;
stimulusTimes = spikes.info.stimtimes{1} + spikes.info.trialonset;
nTrials = length(stimulusTimes); % number of stimulus onsets
if (~isempty(stimulusTimes))
    stimIDs = round(Fs * stimulusTimes);
    if (size(stimIDs,1) > size(stimIDs,2)); stimIDs = stimIDs'; end
    ids = bsxfun(@plus,stimIDs,(sPre + 1 : sPost)');
    ids(ids < 1)    = 1;
    ids(ids > Nmax) = Nmax;
end

for iClus = 1:numclusts

    signalClus = zeros(Nmax, 1);
    id = (spikes.assigns == clusterIDs(iClus));
    spiketimes = spikes.spiketimes(id);
    
    signalClus(round(Fs * spiketimes)) = 1;

    if (stimulusAmp == 0)
        
        N = sum(signalClus) / Tmax; % just return average firing rate
        
    else 
        
        spikeWindows = signalClus(ids); % extract windows around each stimulus onset
        
        N = spikeWindows; % [samples per window x ntrials] 
        N = reshape(N, sBin, nBins, nTrials); % [samples per bin x nbins x ntrials]
                
%         N = sum(N); % sum all spikes in bin [1 x [nbins x ntrials]]
%         N = reshape(N,nBins,[]); % [nbins x ntrials]
%         SpikesBinTrials{clusterIDs(iClus)} = N;
%         
%         N = mean(spikeWindows,2); % sum across trials 
%         N = reshape(N,sBin,[]); % [samples per bin x nbins]
%         N = sum(N); % sum across samples per bin [1 x nbins]
%         fRateTime(:,clusterIDs(iClus)) = N ./ (tBin / 1000);
%         
%         fRate = sum(spikeWindows);
%         fRate = fRate / (tWin / 1000);
%         fRateAmps(clusterIDs(iClus)) = mean(fRate);
    end    
    
    SpikeBinTrials{iClus} = N;
    
end