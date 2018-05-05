function psr_ms_combine_offsets(files,parameters)

delta    = parameters.ms.offset.delta;
nOffsets = length(parameters.ms.offset.min);
nProbes  = size(files,1);

% Initialize
stimOffsetMean = [];
stimOffset     = [];
stimAmps       = [];

for iProbe = 1:nProbes
    metadata = [];
    filename = files{iProbe};
    if (~isempty(filename))
        load(filename,'metadata');
        if (~isfield(metadata,'stimoffsetProbe'))
            if (isfield(metadata,'stimoffset') && isfield(metadata,'stimamps'))
                stimOffset = [stimOffset;metadata.stimoffset];
                stimAmps   = [stimAmps;  metadata.stimamps];
            end
        end
    end
end

I = ~isnan(stimAmps);

if (any(I))
    
    stimOffset = stimOffset(I,:);
    stimAmps   = stimAmps  (I);
    
    n = size(stimOffset,2);
    weights = stimAmps / max(stimAmps);
    weights = repmat(weights,1,n);
    
    stimOffsetAll = stimOffset(:);
    weightsAll    =    weights(:);
    
    dt = parameters.ms.offset.win;
    
    for iOffset = nOffsets:-1:1
        
        t = parameters.ms.offset.min(iOffset);
        
        I = stimOffsetAll < t + dt & stimOffsetAll > t - dt;
        
        stimOffset = stimOffsetAll(I);
        weights    =    weightsAll(I);
        
        if (~isempty(stimOffset))
            
            offset = min(stimOffset);

            subs  = round((stimOffset-offset) / delta) + 1;
            histw = accumarray(subs,weights);
            [~,I] = max(histw);

            stimOffsetMean(iOffset) = delta * (I - 1) + offset;
            
        end
    end
end

% Save output
if (~isempty(stimOffsetMean))
    for iProbe = 1:nProbes
        filename = files{iProbe};
        if (isempty(filename)); continue; end
        load(filename,'metadata');
        metadata.stimoffsetProbe = stimOffsetMean;
        metadata = orderfields(metadata);
        save(filename,'metadata','-append');
    end
end

end