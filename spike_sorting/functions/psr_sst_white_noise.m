function spikes = psr_sst_white_noise(spikes,parameters)

if (isempty_field(spikes,'spikes.clusters.noise.id')); return; end

clustIDs  = [spikes.clusters.noise.id];
nClusts   = length(clustIDs);
nPoints   = size(spikes.waveforms,2);
nChan     = size(spikes.waveforms,3);
precision = 10^parameters.general.precision;
sigma     = spikes.info.bgn;

for iClust = 1:nClusts
    
    spikeIDs = spikes.assigns == clustIDs(iClust);
    chanIdxStruct = spikes.clusters.noise(iClust);
    fields  = fieldnames(chanIdxStruct);
    keep    = ~strcmp(fields,'id');
    fields  = fields(keep);    
    nFields = length(fields);
    I = false(1,nChan);
    
    for iField = 1:nFields
        ids = chanIdxStruct.(fields{iField});
        I(ids) = true;
    end
        
    chanIDs = find(I);
    for iChan = chanIDs
        R = normrnd(0,sigma(iChan),sum(spikeIDs),nPoints);
        spikes.waveforms(spikeIDs,:,iChan) = int16(precision * R);
    end
           
end

end