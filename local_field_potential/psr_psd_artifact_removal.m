function artifacts = psr_lfp_artifact_detection_psd(signal,parameters)

% Artifact removal based on power spectral density

nTrials = length(signal);
fRange  = parameters.lfp.artifact.freqRange;
Fr      = parameters.Fr;
sWin    = parameters.lfp.artifact.window * Fr;
sOff    = 0.5 * sWin; % Take half of window
sSec    = sWin - sOff;

%% Calculate power spectral density for every data section

psdAll = [];
secAll = [];

for iTrial = 1:nTrials

    signalTrial   = mean(signal{iTrial}); % average across channels
    signalSegment = buffer(signalTrial,sWin,sOff);
    signalSegment = signalSegment(:,2:end-1); % ignore first and last segments
    nSecs = size(signalSegment,2);
   
    for iSec = 1:nSecs 
        signalSec = signalSegment(:,iSec);
        psd = pwelch(signalSec,[],[],fRange,Fr);
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
thresh  = parameters.lfp.artifact.threshPSD * psr_mad(psdSum);
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