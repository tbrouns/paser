function [chanIDs,ampAbs,ampRel,p2p,peakLocs] = psr_sst_cluster_amp(spikes, clustID, parameters)

% PSR_SST_CLUSTER_AMP - Calculates amplitude info for cluster
%
% Syntax:  [chanIDs,ampAbs,ampRel,p2p,peakLocs] = psr_sst_cluster_amp(spikes, clustID, parameters)
%
% Inputs:
%    spikes     - See README
%    clustID    - Cluster ID for which we want to calculate the amplitude metrics 
%    parameters - See README and PSR_PARAMETERS_GENERAL
%  
% Outputs:
%    chanIDs  - Channel IDs that cross the spike detection threshold
%    ampAbs   - Absolute amplitude for each channel
%    ampRel   - Relative amplitude for each channel (relative to the spike
%               detection threshold)
%    p2p      - Peak-to-peak amplitude for each channel
%    peakLocs - Location of the peak in each channel

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

threshold  = parameters.cluster.thresh * spikes.info.bgn;
signThresh = sign(parameters.spikes.thresh);
spikeIDs   = ismember(spikes.assigns,clustID);

waveforms = spikes.waveforms(spikeIDs,:,:);
waveforms = psr_int16_to_single(waveforms,parameters);
waveforms = signThresh * squeeze(median(waveforms,1)); % calculate median because of outliers

[ampAbs,peakLocs] = max(waveforms); % maximum amplitude over sample points per channel
chanIDs           = find(ampAbs  > abs(threshold)); % channels that cross threshold with median amplitude
ampRel            = ampAbs ./ abs(threshold);
p2p               = ampAbs - min(waveforms);

end