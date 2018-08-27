function output = psr_wrapper_cSPIKE(spiketrains,parameters)

% spiketrains: A cell array with SpikeTrains{1} containing an
%              array of spike times [spike1 spike2 ...spikeN] for the
%              first spike train and respectively for the other spike
%              trains. The object accepts only spike data aligned as row
%              vectors
%
% trange: [1x2] vector of start and end time of recording
% twin:   [1x2] vector of start and end time of window of interest
% 
% See: http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/cSPIKE.html

rootPath = parameters.analysis.cspike.path;
addpath(rootPath);
addpath([rootPath '\cSPIKEmex']);

trange = parameters.analysis.cspike.trange;
if (~isempty_field(parameters,'parameters.analysis.cspike.twin')); twin = parameters.analysis.cspike.twin;
else,                                                              twin = trange;
end
time1 = twin(1);
time2 = twin(2);

output = [];

% Ignore empty spike trains
keep = ~cellfun(@isempty,spiketrains);
spiketrains = spiketrains(keep);

nTrains = length(spiketrains);

if (nTrains > 1)
    
    % Convert to double
    
    for iTrain = 1:nTrains
        spiketrains{iTrain} = double(spiketrains{iTrain});
    end
    
    % Calculate metrics
    STS = SpikeTrainSet(spiketrains, trange(1), trange(2));
    
    output.ISIdistance                          = STS.ISIdistance                         (time1, time2);
    output.ISIdistanceMatrix                    = STS.ISIdistanceMatrix                   (time1, time2);
    output.ISIdistanceProfile                   = STS.ISIdistanceProfile                  (time1, time2);
    output.AdaptiveISIdistance                  = STS.AdaptiveISIdistance                 (time1, time2);
    output.AdaptiveISIdistanceMatrix            = STS.AdaptiveISIdistanceMatrix           (time1, time2);
    output.AdaptiveISIdistanceProfile           = STS.AdaptiveISIdistanceProfile          (time1, time2);
    output.SPIKEdistance                        = STS.SPIKEdistance                       (time1, time2);
    output.SPIKEdistanceMatrix                  = STS.SPIKEdistanceMatrix                 (time1, time2);
    output.SPIKEdistanceProfile                 = STS.SPIKEdistanceProfile                (time1, time2);
    output.AdaptiveSPIKEdistance                = STS.AdaptiveSPIKEdistance               (time1, time2);
    output.AdaptiveSPIKEdistance                = STS.AdaptiveSPIKEdistanceMatrix         (time1, time2);
    output.AdaptiveSPIKEdistanceProfile         = STS.AdaptiveSPIKEdistanceProfile        (time1, time2);
    output.RateIndependentSPIKEdistance         = STS.RateIndependentSPIKEdistance        (time1, time2);
    output.AdaptiveRateIndependentSPIKEdistance = STS.AdaptiveRateIndependentSPIKEdistance(time1, time2);
    output.SPIKEsynchro                         = STS.SPIKEsynchro                        (time1, time2);
    output.AdaptiveSPIKEsynchro                 = STS.AdaptiveSPIKEsynchro                (time1, time2);
    
    [SPIKESM,SPIKEOM,NormSPIKEOM] = STS.SPIKESynchroMatrix(time1, time2);
    
    output.SPIKESM     = SPIKESM;
    output.SPIKEOM     = SPIKEOM;
    output.normSPIKEOM = NormSPIKEOM;
    
    [SPIKESM,SPIKEOM,NormSPIKEOM] = STS.AdaptiveSPIKESynchroMatrix(time1, time2);
    
    output.AdaptiveSPIKESM     = SPIKESM;
    output.AdaptiveSPIKEOM     = SPIKEOM;
    output.AdaptiveNormSPIKEOM = NormSPIKEOM;
    
    % synchro: a cell array identical to spike train set but instead of spike
    %          times it contains SPIKE-synchronization values of each spike
    % Sorder: a cell array identical to spike train set but instead of spike
    %         times it contains SPIKE-order values of each spike
    % STOrder: a cell array identical to spike train set but instead of spike
    %          times it contains spike-train-order values of each spike
    
    [Synchro,Sorder,STOrder] = STS.SPIKEsynchroProfile(time1, time2);
    
    output.synchro = Synchro;
    output.Sorder  = Sorder;
    output.STOrder = STOrder;
    
    [Synchro,Sorder,STOrder] = STS.AdaptiveSPIKEsynchroProfile(time1, time2);
    
    output.AdaptiveSynchro = Synchro;
    output.AdaptiveSorder  = Sorder;
    output.AdaptiveSTOrder = STOrder;
    
    %% Additional pair-by-pair basis analysis
    
    output.RateIndependentSPIKEdistanceMatrix         = NaN(nTrains,nTrains);
    output.AdaptiveRateIndependentSPIKEdistanceMatrix = NaN(nTrains,nTrains);
    
    for iTrain = 1:nTrains
        for jTrain = 1:nTrains
            spiketrain_pair = [spiketrains(iTrain) spiketrains(jTrain)];
            STS_pair = SpikeTrainSet(spiketrain_pair, trange(1), trange(2));
            output.RateIndependentSPIKEdistanceMatrix        (iTrain,jTrain) = STS_pair.RateIndependentSPIKEdistance        (time1, time2);
            output.AdaptiveRateIndependentSPIKEdistanceMatrix(iTrain,jTrain) = STS_pair.AdaptiveRateIndependentSPIKEdistance(time1, time2);
        end
    end
    
    %% Output
    output.trainIDs = find(keep)';
    output = orderfields(output);
    
end

end