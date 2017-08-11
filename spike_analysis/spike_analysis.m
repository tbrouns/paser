function spike_analysis(spikes,clusters)

clusterIDs = cell2mat({clusters.vars.id});
clustFlags = cell2mat({clusters.vars.flag});
clusterIDs(clustFlags) = [];

numclusts = length(clusterIDs);

for iClust = 1:numclusts
    id = spikes.assigns == clusterIDs(iClust);
    spikes = ss_spike_removal(spikes,id);
end

clusterIDs  = cell2mat({clusters.vars.id});
clusterIDs(~clustFlags) = [];
numclusts = length(clusterIDs);

stimulusTimes = spikes.artifacts;
nstimuli = length(stimulusTimes);

Tmax = max(spikes.spiketimes);
Fs   = spikes.params.Fs;
Fr   = 25; % Hz

twin = 250; 
time_window = Fs*twin; % ms

signalStim = zeros(Fs*Tmax,1);
signalStim(round(Fs*stimulusTimes)) = 1;

figure
for iClus = 1:numclusts   
    clf
    signalClus = zeros(Fs*Tmax,1);
    id = spikes.assigns == clusterIDs(iClus);
    spiketimes = spikes.spiketimes(id);
    signalClus(round(Fs*spiketimes)) = 1;
    xc = xcorr(signalStim,signalClus,time_window,'coeff');
    
    N = mod(length(xc),(Fs / Fr));
    xc = xc(1:end-N);
    xc = reshape(xc,[],(Fs / Fr));
    xc = mean(xc,2);
    
    t = linspace(-twin,twin,length(xc));
    
    bar(t,xc,1.0);
end

% for iStim = 1:nstimuli
%     id = find(spikes.spiketimes > stimulusTimes(iStim) && spikes.spiketimes <= stimulusTimes(iStim) + time_window / 1000);
%     spiketimes = spikes.spiketimes(id);
%     assigns    = spikes.assigns(id);
%     
% end
% 
% for iclust = 1:numclusts
%     
%     spike_times = spikes.spiketimes(spikes.assigns == clusterIDs(iclust));
%         
% end

end