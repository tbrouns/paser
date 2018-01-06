function correlations = psr_sst_cluster_corr(spikes,parameters,filterSpikes)

if (nargin < 3); filterSpikes = true; end

nClust = max(spikes.assigns);

% Filter spikes
if (isfield(spikes,'delete') && filterSpikes)
    spikes = psr_sst_filter_spikes(spikes,parameters,'delete');
end

% Convert to single
if (isa(spikes.waveforms,'int16')); spikes.waveforms = psr_single(spikes.waveforms,parameters); end

% Calculate features of clusters
correlations = NaN(nClust,nClust);

for iClust = 1:nClust
    
    % Extract cluster ID
    nspikes = sum(spikes.assigns == iClust);
    if (nspikes == 0); continue; end
    
    % Pair-wise cluster correlations
    for jClust = (iClust+1):nClust
        correlations(iClust,jClust) = clusterCorr(spikes,iClust,jClust,parameters);
    end
end

end

function r = clusterCorr(spikes, clustID_1, clustID_2, parameters)

maxlag = round(0.001 * parameters.spikes.max_desync * spikes.Fs);

spikeIDs_1 = ismember(spikes.assigns, clustID_1);
spikeIDs_2 = ismember(spikes.assigns, clustID_2);

if (isempty(spikeIDs_2)); r = NaN; return; end

waves_1 = spikes.waveforms(spikeIDs_1,:);
waves_2 = spikes.waveforms(spikeIDs_2,:);

waves_1 = squeeze(mean(waves_1,1));
waves_2 = squeeze(mean(waves_2,1));

r = xcorr(waves_1,waves_2,maxlag,'coeff');
r = max(r);

end