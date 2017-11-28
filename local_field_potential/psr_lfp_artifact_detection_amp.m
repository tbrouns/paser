function artifacts = psr_lfp_artifact_detection_amp(signal,parameters)

% Artifact removal based on z-score amplitude

nTrials = length(signal);
Fs = parameters.Fr; 
Fr = parameters.lfp.artifact.freq;
[~,slength] = cellfun(@size,signal);
slength = [0;cumsum(slength)];

signal = cat(2,signal{:}); % Combine data across all stimulus conditions
zscore = bsxfun(@minus,  signal,mean(signal,2)); % Calculate z-score
zscore = bsxfun(@rdivide,zscore, std(zscore,[],2));
zscore = abs(mean(zscore)); % Average across channels

tsection  = parameters.lfp.artifact.tsection;
threshold = parameters.lfp.artifact.threshAmp;
acutoff   = findThreshold(zscore,threshold,tsection,Fs);

t     = (1:size(zscore,2)) / Fs;
zfilt = gaussfilt(t,zscore,1/Fr);

[pks_max,locs_max_all] = findpeaks(double(zscore)); % detect peaks in raw data signal
[      ~,locs_min_all] = findpeaks(double(-zfilt)); % detect on- and offset points in filtered data

% Find peaks above threshold

locs_max = locs_max_all(pks_max > acutoff);

% Find adjacent minima

nsamples  = size(signal,2);
artifacts = findRange(locs_max,locs_min_all,nsamples);
artifacts = artifacts';
artifactsAll = cell(nTrials,1);

for iTrial = 1:nTrials
    smin = slength(iTrial);
    smax = slength(iTrial + 1);
    id = artifacts(:,2) > smin & artifacts(:,1) < smax;
    artifactsTrial = artifacts(id,:) - smin;
    artifactsTrial(artifactsTrial < 1)    = 1;
    artifactsTrial(artifactsTrial > smax) = smax;
    artifactsTrial = unique(artifactsTrial,'rows');
    artifactsAll{iTrial} = (artifactsTrial - 1) / Fs;
end

artifacts = artifactsAll;

end

function acutoff = findThreshold(data,threshold,tsection,fs)

% Find background noise

nlength  = size(data,2);
nsamples = round(tsection * fs); % cut data in sections of X seconds
nstep    = round(0.1 * (nsamples)); % move window with steps of 0.1*nsection
stdev    = [];

iStart = 1;
iEnd   = iStart + nsamples;
STOP   = 0;

while (~STOP)
    
    if (iEnd > nlength)
        iStart = nlength - nsamples; 
        iEnd   = nlength;
        STOP   = 1;
    end
    
    if (  iEnd > nlength); iEnd   = nlength; end
    if (iStart <       1); iStart = 1;       end
    
    data_section = data(:,iStart:iEnd);
    stdev  = [stdev;median(abs(data_section)) / 0.6745]; %#ok
    iStart = iStart + nstep;
    iEnd   = iStart + nsamples;
end
    
stdev   = min(stdev);
acutoff = threshold * stdev;

end

function locs_artifact = findRange(locs,locs_all,nsamples)

N             = length(locs);
locs_artifact = zeros(2,N);

for iLoc = 1:N
    loc_max = locs(iLoc);
    id = (loc_max - locs_all) > 0;
    [~,I1] = min(loc_max - locs_all( id));
    [~,I2] = max(loc_max - locs_all(~id));
    n = find(~id,1) - 1;
    if (~isempty(I1)); locs_artifact(1,iLoc) = locs_all(I1);
    else,              locs_artifact(1,iLoc) = 1;
    end
    if (~isempty(I2)); locs_artifact(2,iLoc) = locs_all(n + I2);
    else,              locs_artifact(2,iLoc) = nsamples;
    end
end

end