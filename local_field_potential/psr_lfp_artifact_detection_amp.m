function artifacts = psr_lfp_artifact_detection_amp(data,parameters)

[data,nBlocks] = psr_lfp_conversion(data);

% Set parameters

Fs                = parameters.Fr;
tSection          = parameters.lfp.artifact.amp.tsection;
threshFactorUpper = parameters.lfp.artifact.amp.upper;
threshFactorLower = parameters.lfp.artifact.amp.lower;
sWinSlope         = round(Fs * parameters.lfp.artifact.amp.slope / 1000);

[~,sLength] = cellfun(@size,data);
sLength     = [0,cumsum(sLength)];

% Combine data across all blocks

if (parameters.lfp.artifact.amp.cat); data = {cell2mat(data)}; end

% Find artifacts

artifacts = [];
mBlocks = length(data);

for iBlock = 1:mBlocks
        
    % Amplitude 
    
    dataBlock = data{iBlock};
    dataBlock = nanmean(abs(dataBlock),1); % Average across channels
    
    artifactsTrial = findArtifacts(dataBlock,tSection,Fs,threshFactorUpper,threshFactorLower);
    artifacts = [artifacts;artifactsTrial+sLength(iBlock)];
    
    % Derivative
    
    derivative = abs(dataBlock(1+sWinSlope:end) - dataBlock(1:end-sWinSlope));
    
    artifactsTrial = findArtifacts(derivative,tSection,Fs,threshFactorUpper,threshFactorLower);
    artifacts = [artifacts;artifactsTrial+sLength(iBlock)];
    
end

% Find artifacts for each trial

artifactsAll = cell(nBlocks,1);

for iBlock = 1:nBlocks
    smin = sLength(iBlock);
    smax = sLength(iBlock + 1);
    id = artifacts(:,2) > smin & artifacts(:,1) < smax;
    artifactsTrial = artifacts(id,:) - smin;
    artifactsTrial(artifactsTrial < 1)    = 1;
    artifactsTrial(artifactsTrial > smax) = smax;
    artifactsTrial = unique(artifactsTrial,'rows');
    artifactsAll{iBlock} = (artifactsTrial - 1) / Fs;
end

artifacts = artifactsAll;

end

function bgNoise = findBackgroundNoise(data,tSection,Fs)

% Find background noise

% We determine the median absolute deviation (MAD) in data sections of length
% "tSection" and then use the smallest MAD as the 

nLength  = length(data);
nSamples = round(tSection * Fs); % cut data in sections of X seconds
nStep    = round(0.1 * (nSamples)); % move window with steps of 0.1*nsection
bgNoise  = [];

iStart = 1;
iEnd   = iStart + nSamples;
STOP   = 0;

while (~STOP)
    
    if (iEnd > nLength)
        iStart = nLength - nSamples;
        iEnd   = nLength;
        STOP   = 1;
    end
    
    if (  iEnd > nLength); iEnd   = nLength; end
    if (iStart <       1); iStart = 1;       end
    
    data_section = data(iStart:iEnd);
    bgNoise  = [bgNoise;psr_mad(data_section)]; %#ok
    iStart = iStart + nStep;
    iEnd   = iStart + nSamples;
end

bgNoise = min(bgNoise);

end

function artifacts = findArtifacts(signal,tSection,Fs,threshFactorUpper,threshFactorLower)

sLength = length(signal);
bgNoise = findBackgroundNoise(signal,tSection,Fs);
threshUpper = threshFactorUpper * bgNoise;
threshLower = threshFactorLower * bgNoise;

[peakAmps,peakLocs] = findpeaks(double(signal)); % detect peaks in raw data signal
peakLocs = peakLocs(peakAmps > threshUpper); % Find peaks above threshold
ids = double(signal < threshLower);
ids(peakLocs) = 2;
ids(1)        = 1;
ids(sLength)  = 1;
IDs = find(ids);

locs    = ids(ids > 0);
peakIDs = find(locs == 2);
onsets  = IDs(peakIDs - 1);
offsets = IDs(peakIDs + 1);
onsets  = onsets (ids(onsets)  == 1);
offsets = offsets(ids(offsets) == 1);
artifacts = ([onsets;offsets])';

end