function [fRateStimTime,fRateStimAmps] = createPSTH(spikes,clusters,clusterMax,tPre,tPost,tbin)

twin  = tPre + tPost;
Tmax  = spikes.info.detect.dur;
Fs    = spikes.params.Fs;
sPre  = Fs*(tPre /1000);
sPost = Fs*(tPost/1000);
sbin  = Fs*(tbin /1000);

if (~isfield(spikes,'artifacts')); control = 1;
else                               control = 0;
end

clusterIDs = cell2mat({clusters.vars.id});
numclusts  = length(clusterIDs);

if (~control);
    stimulusTimes = spikes.artifacts;
    stimIDs = round(Fs*stimulusTimes);
    if (size(stimIDs,1) > size(stimIDs,2)); stimIDs = stimIDs'; end
    ids = bsxfun(@plus,stimIDs,(-sPre+1:sPost)');
    ids(ids < 1)              = 1;
    ids(ids > floor(Fs*Tmax)) = floor(Fs*Tmax);
end
    
fRateStimTime = zeros(twin/tbin,clusterMax);
fRateStimAmps = NaN(1,clusterMax);

for iClus = 1:numclusts

    signalClus = zeros(floor(Fs*Tmax),1);
    id = (spikes.assigns == clusterIDs(iClus));
    spiketimes = spikes.spiketimes(id);
    
%     if (~control); spiketimes = stimulusTimes + (100 * rand(size(stimulusTimes)) + 50) / 1000; end
    signalClus(round(Fs*spiketimes)) = 1;

    if (control)
        fRateAvg = sum(signalClus) / Tmax;
        fRateStimTime(:,clusterIDs(iClus)) = fRateAvg;
        fRateStimAmps(  clusterIDs(iClus)) = fRateAvg;
    else
        fRate = mean(signalClus(ids),2);
        fRate = reshape(fRate,sbin,[]);
        fRate = sum(fRate);
        fRateStimTime(:,clusterIDs(iClus)) = fRate ./ (tbin / 1000);

        fRate = sum(signalClus(ids));
        fRate = fRate / (twin / 1000);
        fRateStimAmps(clusterIDs(iClus)) = mean(fRate);
    end        
end