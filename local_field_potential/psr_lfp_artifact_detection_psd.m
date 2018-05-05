function artifacts = psr_lfp_artifact_detection_psd(data,parameters)

[data,nTrials] = psr_lfp_conversion(data);

% Artifact removal based on power spectral density

fRange  = parameters.lfp.artifact.psd.frange;
Fr      = parameters.Fr;
sWin    = parameters.lfp.artifact.psd.win * Fr;
sOff    = 0.5 * sWin; % Take half of window
sSec    = sWin - sOff;

%% Calculate power spectral density for every data section

psdAll = [];
secAll = [];

for iTrial = 1:nTrials

    dataTrial   = data{iTrial};
    if (isfield(dataTrial,'trial')); dataTrial = dataTrial.trial{1}; end
    dataTrial   = mean(dataTrial); % average across channels
    dataSegment = buffer(dataTrial,sWin,sOff);
    dataSegment = dataSegment(:,2:end-1); % ignore first and last segments
    nSecs = size(dataSegment,2);
   
    for iSec = 1:nSecs 
        dataSec = dataSegment(:,iSec);
        psd = pwelch(dataSec,[],[],fRange,Fr);
        if (iSec == 1); psdTrial = zeros(length(psd),nSecs,'single'); end
        psdTrial(:,iSec) = single(psd);
    end
   
    offsets  = sSec * (1:nSecs);
    onsets   = offsets - sSec + 1;
    secTrial = [onsets;offsets;iTrial*ones(size(onsets))];
        
    psdAll = [psdAll,psdTrial]; %#ok
    secAll = [secAll,secTrial]; %#ok

end

psdMean = mean(psdAll,2);
psdAll  = bsxfun(@rdivide,psdAll,psdMean); % Normalize

% Total PSD threshold
psdSum  = sum(psdAll); 
thresh  = parameters.lfp.artifact.psd.thresh * psr_mad(psdSum);
artifacts = psdSum > thresh;

% Extract intervals
artifacts = secAll(:,artifacts);
artifactsAll = cell(nTrials,1);
for iTrial = 1:nTrials
    artifactsTrial = artifacts(:,artifacts(3,:) == iTrial);
    timings = artifactsTrial(1,:);
    if (~isempty(timings)) 
        d  = diff(timings); % Combine adjacent intervals
        offsets = find(d > sOff); 
        onsets  = offsets + 1;
        offsets = artifactsTrial(2,unique([offsets length(timings)]));
        onsets  = artifactsTrial(1,unique([1 onsets]));
        artifactsAll{iTrial} = (([onsets;offsets] - 1) / Fr)';
    end
end

artifacts = artifactsAll; % output

end