function artifacts = psr_lfp_artifact_detection_psd(data,parameters)

[data,nBlocks] = psr_lfp_conversion(data);

% Artifact removal based on power spectral density

fRange  = parameters.lfp.artifact.psd.frange;
Fr      = parameters.Fr;
sWin    = parameters.lfp.artifact.psd.win * Fr;
sOff    = 0.5 * sWin; % Take half of window
sSec    = sWin - sOff;

%% Calculate power spectral density for every data section

psdAll = cell(0,0);
secAll = cell(0,0);

for iBlock = nBlocks:-1:1

    dataBlock = data{iBlock};
    if (isfield(dataBlock,'trial')); dataBlock = dataBlock.trial{1}; end
    
    nChans = size(dataBlock,1);
    if (~isempty(dataBlock))
        
        for iChan = 1:nChans
                        
            dataChan    = dataBlock(iChan,:); % average across channels
            dataSegment = buffer(dataChan,sWin,sOff); % extract the segments, with overlap
            dataSegment = dataSegment(:,2:end-1); % ignore first and last segments
            nSecs = size(dataSegment,2); 

            psdTemp = [];
            for iSec = nSecs:-1:1 
                dataSec = dataSegment(:,iSec);
                psd = pwelch(dataSec,[],[],fRange,Fr);
                if (isempty(psdTemp)); psdTemp = zeros(length(psd),nSecs,'single'); end % Initialize
                psdTemp(:,iSec) = single(psd);
            end

            offsets = sSec * (1:nSecs);
            onsets  = offsets - sSec + 1;
            secTemp = [onsets;offsets;iBlock*ones(size(onsets))];

            % Save
            psdAll{iChan,iBlock} = psdTemp;
            secAll{iChan,iBlock} = secTemp;
        end
    end
end

keepPsd = any(~cellfun(@isempty,psdAll),1);
keepSec = any(~cellfun(@isempty,secAll),1);

psdAll = psdAll(:,keepPsd);
secAll = secAll(:,keepSec);

psdAll = cell2mat(psdAll);
secAll = cell2mat(secAll);

psdMean = mean(psdAll,2); % Average across all segments
psdAll  = bsxfun(@rdivide,psdAll,psdMean); % Normalize

% Total PSD threshold
psdSum  = sum(psdAll); % Sum across all frequencies 
thresh  = parameters.lfp.artifact.psd.thresh * psr_mad(psdSum); % Artifact thresholding
artifacts = psdSum > thresh;

% Extract intervals
artifacts = secAll(:,artifacts);
artifactsAll = cell(nBlocks,1);
for iBlock = 1:nBlocks
    artifactsBlock = artifacts(:,artifacts(3,:) == iBlock);
    timings = artifactsBlock(1,:);
    if (~isempty(timings)) 
        d = diff(timings); % Combine adjacent intervals
        offsets = find(d > sOff); 
        onsets  = offsets + 1;
        offsets = artifactsBlock(2,unique([offsets length(timings)]));
        onsets  = artifactsBlock(1,unique([1 onsets]));
        artifactsAll{iBlock} = (([onsets;offsets] - 1) / Fr)';
    end
end

artifacts = artifactsAll; % output

end