function MSR = psr_multiscale_relevance(spiketimes,tobs,parameters)

% PSR_MULTISCALE_RELEVANCE - Calculates the multi-scale relevance metric
% 
% References: 
% [1] "Finding informative neurons in the brain using Multi-Scale
%     Relevance", Cubero et al. (2018)
% 
% Syntax:  MSR = psr_multiscale_relevance(spiketimes,tobs,parameters)
%
% Inputs:
%    spiketimes - Spike times [in sec] for a particular neuron, given as a vector
%    tobs       - Total observation time [sec]
%    parameters - See PSR_PARAMETERS_ANALYSIS
%
% Outputs:
%    MSR - The multi-scale relevance measure for the spike train

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

dt     = parameters.analysis.msr.dt_min; % Minimum time scale 
dt_min = log10(dt);
dt_max = log10(tobs);
dTs    = logspace(dt_min,dt_max,parameters.analysis.msr.dt_step); % Time scale range [sec]

% Initialize
Ndt = length(dTs);
Hs  = NaN(Ndt,1);
HK  = NaN(Ndt,1);
MSR = NaN;

Fs = 10 / dt; % Make sampling frequency 10 times greater than bin frequency
spikepoints = round(spiketimes * Fs);
tobs        = round(      tobs * Fs);
edges = 0.5:tobs+0.5;
spikevector = histcounts(spikepoints,edges);
M = sum(spikevector);
if (M == 0); return; end

for idt = 1:Ndt
       
    dt = dTs(idt);
    
    % Resolution
    v = ones(1,round(dt * Fs));
    ks = conv(spikevector,v,'same');
    n = length(v);
    phase = round(0.5 * n);
    ksd = ks;
    if (phase < n);   ksd = downsample(ks,n,phase); end
    if all(ksd == 0); ksd = downsample(ks,n,n-1);   end
    Hs(idt) = -nansum((ksd/M).*log(ksd/M)/log(M));
    
    % Relevance
    edges = 0.5:max(ksd)+0.5;
    if (length(edges) >= 2)
        mk = histcounts(ksd,edges);
        k  = edges(1:end-1) + 0.5;
        HK(idt) = -nansum((k.*mk/M).*log(k.*mk/M)/log(M));
    end
end

del = isnan(HK);
HK(del) = [];
Hs(del) = [];

% Area under curve (MSR)
[Hs,Isort] = sort(Hs);
HK = HK(Isort);
MSR = trapz(Hs,HK);

end