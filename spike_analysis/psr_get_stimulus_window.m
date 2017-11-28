function [SpikeBinTrials,SpikeTimesTrial] = psr_get_stimulus_window(spikes,params)

tPre  = params.t_win(1);
tPost = params.t_win(end);
tBin  = params.t_bin;
tDel  = params.t_del;

tWin  = tPost - tPre;
Tmax  = spikes.info.dur;
Fs    = spikes.Fs;
Nmax  = floor(Fs * Tmax) + 1;
nBins = tWin / tBin;
sPre  = Fs * (tPre  / 1000); % pre- stimulus window
sPost = Fs * (tPost / 1000); % post-stimulus window
sBin  = Fs * (tBin  / 1000); % bin size in stimulus window
sDel  = Fs * (tDel  / 1000); % spike deletion window

clusterIDs      = [spikes.clusters.vars.id];
nClusts         = length(clusterIDs);
SpikeBinTrials  =  cell(1,nClusts);
SpikeTimesTrial =  cell(2,nClusts);

% Extract stimulus windows

stimulusTimes = spikes.info.stimtimes{1}(:,1) + (params.stimOffset / 1000);
nTrials = length(stimulusTimes); % number of stimulus onsets

if (~isempty(stimulusTimes))
    
    stimulusTimes = round(Fs * stimulusTimes) + 1;
    if (size(stimulusTimes,1) > size(stimulusTimes,2)); stimulusTimes = stimulusTimes'; end
    
    del = bsxfun(@plus,stimulusTimes,(sDel(1) : sDel(2) - 1)');
    del(del < 1)    = 1;
    del(del > Nmax) = Nmax;
    del = del(:);
    
    ids = bsxfun(@plus,stimulusTimes,(sPre + 1 : sPost)');
    ids(ids < 1)    = 1;
    ids(ids > Nmax) = Nmax;
        
    d = round(Fs * (params.acf_win / 1000));
    preId = bsxfun(@plus,stimulusTimes,(-d : 0)');
    pstId = bsxfun(@plus,stimulusTimes,( 0 : d)');
    preId = preId(:);
    pstId = pstId(:);
    preId(preId < 1)    = [];
    pstId(pstId < 1)    = [];
    preId(preId > Nmax) = [];
    pstId(pstId > Nmax) = [];

    for iClus = 1:nClusts

        % Binary spiking vector
        signal = false(Nmax, 1);
        id = (spikes.assigns == clusterIDs(iClus));
        spiketimes = spikes.spiketimes(id) - spikes.info.trialonset;
        spiketimes = round(Fs * spiketimes) + 1; % in sample number
        signal(spiketimes) = true;
        signal(del)        = false;
        
        SpikeBinCount = signal(ids); % extract windows around each stimulus onset [samples per window x ntrials]
        SpikeBinCount = reshape(SpikeBinCount, sBin, nBins, nTrials); % [samples per bin x nbins x ntrials]        
        SpikeBinTrials{iClus} = SpikeBinCount;

        % Pre- and post-stimulus spiking vectors

        signalPre = false(Nmax,1);
        signalPst = false(Nmax,1);

        SpikeTimesPre = signal;
        SpikeTimesPst = signal;

        signalPre(preId) = true;
        signalPst(pstId) = true;

        SpikeTimesPre(~signalPre) = false;
        SpikeTimesPst(~signalPst) = false;

        SpikeTimesPre = (find(SpikeTimesPre) - 1) / Fs; 
        SpikeTimesPst = (find(SpikeTimesPst) - 1) / Fs; 

        SpikeTimesTrial{1,iClus} = SpikeTimesPre;
        SpikeTimesTrial{2,iClus} = SpikeTimesPst;
        SpikeTimesTrial{3,iClus} = Tmax;

    end

end

end