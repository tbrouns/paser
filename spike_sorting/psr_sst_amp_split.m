function [amplitudes,spikeIDs] = psr_sst_amp_split(spikes,clustID,parameters)

binSpikes      = parameters.cluster.split.bin;
binSmooth      = parameters.cluster.split.smooth;
threshProm     = parameters.cluster.split.prom;
threshWidth    = parameters.cluster.split.width;
threshTerminal = parameters.cluster.split.terminal;

thresh   = abs(mean(spikes.info.thresh));
minBin   = 10^(-parameters.general.precision + 1) / thresh; % limit bin size by precision of data
spikeIDs = find(spikes.assigns == clustID);
ampClust = psr_sst_amp(spikes,clustID,parameters);
amplitudes = ampClust; % Initialize
clear waves spikes;

if (~isempty(amplitudes))
      
    nspikes = length(ampClust);
    bin     = binSpikes * range(ampClust) / nspikes;
    if (bin < minBin); bin = minBin; end
    edges   = min(ampClust):bin:max(ampClust);
    
    if (length(edges) >= 3) % condition needed for findpeaks
        
        counts = histcounts(ampClust,edges); % amplitude distribution
        counts = psr_gauss_smoothing(counts,binSmooth);
        counts = counts / max(counts); % normalize
        nbins  = length(counts);
        counts = [0,counts]; % to ensure that first bin can also be peak
        
        [pks,locs,widths,proms] = findpeaks(counts);
        
        % Search left and right of maximum amplitude
        
        [~,I]   = max(pks);
        locMax  = locs(I);
        
        idLeft = locs < locMax;
        idRght = locs > locMax;
        
        locsLeft  =   locs(idLeft);
        locsRght  =   locs(idRght);
        promsLeft =  proms(idLeft);
        promsRght =  proms(idRght);
        widthLeft = widths(idLeft);
        widthRght = widths(idRght);
        
        pkIdLeft = find(promsLeft > threshProm & widthLeft > threshWidth,1,'last');
        pkIdRght = find(promsRght > threshProm & widthRght > threshWidth,1,'first');
        
        locMinLeft = [];
        locMinRght = [];
        
        if (~isempty(pkIdLeft)); locMinLeft = findMinimum(counts,locsLeft(pkIdLeft),locMax); end
        if (~isempty(pkIdRght)); locMinRght = findMinimum(counts,locMax,locsRght(pkIdRght)); end
        
        countsLeft = counts(1:locMax-1);
        countsRght = counts(locMax+1:end);
        
        locTerminalLeft = find(countsLeft < threshTerminal,1,'last');
        locTerminalRght = find(countsRght < threshTerminal,1,'first') + locMax;
        
        locLeft = max([locMinLeft locTerminalLeft]);
        locRght = min([locMinRght locTerminalRght]);
        
        if (isempty(locLeft) || locLeft < 2); locLeft = 2;         end
        if (isempty(locRght));                locRght = nbins + 1; end
        
        ampLeft = edges(locLeft - 1); % Subtract 1 because we added zero to 'counts'
        ampRght = edges(locRght);     % Take next edge
        
        spikeIDs_sub = find(ampClust >= ampLeft & ampClust <= ampRght); % Sub-cluster spike IDs
        amplitudes = ampClust(spikeIDs_sub);
        spikeIDs   = spikeIDs(spikeIDs_sub);
        
    end
    
end

end

function x = psr_gauss_smoothing(x,n)

g = gausswin(n); % create Gaussian smoothing window
g = g / sum(g); % normalize
x = conv(x, g, 'same');

end

function locMin = findMinimum(counts,loc_1,loc_2)

[~,I] = min(counts(loc_1:loc_2)); % Minimum in-between peaks
locMin = I + loc_1 - 1;

end