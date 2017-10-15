function [SpikeBinTrials,FiringRate] = psr_analysis_clusterdata(spikes,params)

tPre  = params.t_array(1);
tPost = params.t_array(end);
tBin  = params.t_bin;
tDel  = params.t_del;

tWin  = tPost - tPre;
Tmax  = spikes.info.detect.dur;
Fs    = spikes.params.Fs;
Nmax  = floor(Fs * Tmax);
nBins = tWin / tBin;
sPre  = Fs * (tPre  / 1000); % pre- stimulus window
sPost = Fs * (tPost / 1000); % post-stimulus window
sBin  = Fs * (tBin  / 1000); % bin size in stimulus window
sDel  = Fs * (tDel  / 1000); % spike deletion window

clusterIDs = [spikes.clusters.vars.id];
numclusts  = length(clusterIDs);
SpikeBinTrials =  cell(1,numclusts);
FiringRate     = zeros(1,numclusts);

% Extract stimulus windows

stimulusAmp   = spikes.info.stimulus;
stimulusTimes = spikes.info.stimtimes{1};
nTrials = length(stimulusTimes); % number of stimulus onsets
if (~isempty(stimulusTimes))
    stimulusTimes = round(Fs * stimulusTimes);
    if (size(stimulusTimes,1) > size(stimulusTimes,2)); stimulusTimes = stimulusTimes'; end
    
    del = bsxfun(@plus,stimulusTimes,(-sDel : sDel)');
    del(del < 1)    = 1;
    del(del > Nmax) = Nmax;
    del = del(:);
    
    ids = bsxfun(@plus,stimulusTimes,(sPre + 1 : sPost)');
    ids(ids < 1)    = 1;
    ids(ids > Nmax) = Nmax;
end

for iClus = 1:numclusts
    
    signalClus = zeros(Nmax, 1);
    id = (spikes.assigns == clusterIDs(iClus));
    spiketimes = spikes.spiketimes(id) - spikes.info.trialonset;
    signalClus(round(Fs * spiketimes)) = 1;
    signalClus(del)                    = 0;
    
    FiringRate(iClus) = sum(signalClus) / Tmax; % just return average firing rate [spikes / sec]
    
    if (stimulusAmp ~= 0)
        
        N = signalClus(ids); % extract windows around each stimulus onset [samples per window x ntrials]
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
        
        SpikeBinTrials{iClus} = N;
        
    end
    
end

end