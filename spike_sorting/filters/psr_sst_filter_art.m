function del = psr_sst_filter_art(spikes)

% Spike removal parameters

Tmax     = sum(spikes.info.dur);
Fs       = spikes.Fs;
nBlocks  = size(spikes.info.artifacts.raw,1);
trialonsets = [0;cumsum(spikes.info.dur)];

% Remove spikes that correspond with artefact locations

artifactsAll = [];
for iBlock = 1:nBlocks
    artifactsAll{iBlock} = spikes.info.artifacts.raw{iBlock,1} + trialonsets(iBlock);
end
if (~isempty(artifactsAll)); artifactsAll = cat(1, artifactsAll{:}); end
artifactsAll = round(Fs * artifactsAll) + 1;
artifactsAll = unique(artifactsAll,'rows');
nArtifacts   =   size(artifactsAll,1);
artifactIDs = [];
for iArtifact = 1:nArtifacts
   I = artifactsAll(iArtifact,1) : artifactsAll(iArtifact,2);
   artifactIDs = [artifactIDs;I'];
end
artifactIDs = unique(artifactIDs);

spikePoints = round(Fs * spikes.spiketimes) + 1; % Get spike points

del = psr_get_spike_ids(spikePoints,artifactIDs); % Spikes inside stimulus window

end