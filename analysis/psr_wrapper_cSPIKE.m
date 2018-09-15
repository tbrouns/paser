function output = psr_wrapper_cSPIKE(spiketrains,parameters)

% PSR_WRAPPER_CSPIKE - Wrapper function for the cSPIKE toolbox
% Calculates a number of different spike train distances 
% 
% References:
% [1] http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/cSPIKE.html
% 
% Syntax:  output = psr_wrapper_cSPIKE(spiketrains,parameters)
%
% Inputs:
%    spiketrains - A cell array with each cell containing an
%                  array of spike times [spike1 spike2 ...spikeN]. 
%                  The object accepts only spike data aligned as row vectors
% 
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    output - Structure of containing many different spike distance metrics
%
% Dependencies: Requires cSPIKE toolbox available from Ref. [1]
%
% See also: SPIKETRAINSET

% These codes are free of charge for research and education purposes
% only. Any commercial or military use of this software is prohibited.
% 
% The software on this site is provided "as-is," without any expressed or
% implied warranty. In no event am I or my host institution liable for any
% damages arising from the use of the software. Since it is distributed for
% free, I do also not take responsibility for any eventual error in it.
% 
% BSD license:
% 
% Copyright (c) 2016 Eero Satuvuori All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% * Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% * Neither the name of the author nor the names of its contributors may be
% used to endorse or promote products derived from this software without
% specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% spiketrains: 
%
% trange: [1x2] vector of start and end time of recording
% twin:   [1x2] vector of start and end time of window of interest

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