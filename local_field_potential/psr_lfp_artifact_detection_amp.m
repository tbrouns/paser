function artifacts = psr_lfp_artifact_detection_amp(signal,parameters)

% Set parameters

Fs = parameters.Fr;
nTrials = length(signal);
tSection = parameters.lfp.artifact.tsection;
threshFactorUpper = parameters.lfp.artifact.threshAmpUpper;
threshFactorLower = parameters.lfp.artifact.threshAmpLower;
sWinSlope  = round(Fs * parameters.lfp.artifact.tSlope / 1000);

[~,slength] = cellfun(@size,signal);
slength = [0;cumsum(slength)];

if (parameters.lfp.artifact.cat) % Combine data across all stimulus conditions
    signal = {cat(2,signal{:})};
end

artifacts = [];

for iTrial = 1:length(signal)
        
    %% Amplitude 
    
    signalTrial = signal{iTrial};
    signalTrial = abs(mean(signalTrial)); % Average across channels
    
    artifactsTrial = findArtifacts(signalTrial,tSection,Fs,threshFactorUpper,threshFactorLower);
    artifacts = [artifacts;artifactsTrial+slength(iTrial)];
    
    %% Derivative
    
    derivative = abs(signalTrial(1+sWinSlope:end) - signalTrial(1:end-sWinSlope));
    
    artifactsTrial = findArtifacts(derivative,tSection,Fs,threshFactorUpper,threshFactorLower);
    artifacts = [artifacts;artifactsTrial+slength(iTrial)];
    
end

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

function stdev = findStd(data,tSection,Fs)

% Find background noise

nLength  = length(data);
nSamples = round(tSection * Fs); % cut data in sections of X seconds
nStep    = round(0.1 * (nSamples)); % move window with steps of 0.1*nsection
stdev    = [];

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
    stdev  = [stdev;psr_mad(data_section)]; %#ok
    iStart = iStart + nStep;
    iEnd   = iStart + nSamples;
end

stdev = min(stdev);

end

% function locs_artifact = findRange(locPeaks,locBounds,nSamples)
% 
% N             = length(locPeaks);
% locs_artifact = zeros(2,N);
% 
% for iLoc = 1:N
%     loc_max = locPeaks(iLoc);
%     id = (loc_max - locBounds) > 0;
%     [~,I1] = min(loc_max - locBounds( id));
%     [~,I2] = max(loc_max - locBounds(~id));
%     n = find(~id,1) - 1;
%     if (~isempty(I1)); locs_artifact(1,iLoc) = locBounds(I1);
%     else,              locs_artifact(1,iLoc) = 1;
%     end
%     if (~isempty(I2)); locs_artifact(2,iLoc) = locBounds(n + I2);
%     else,              locs_artifact(2,iLoc) = nSamples;
%     end
% end
% 
% end

function artifacts = findArtifacts(signal,tSection,Fs,threshFactorUpper,threshFactorLower)

sLength = length(signal);
stdev = findStd(signal,tSection,Fs);
threshUpper = threshFactorUpper * stdev;
threshLower = threshFactorLower * stdev;

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

%% Visualization
% figure;
% artifactsTrial = unique(artifacts,'rows');
% t = ((1:size(signal,2)) - 1) / Fs;
% plot(t,signal); hold on;
% scatter((artifactsTrial(:,1)-1)/Fs,threshLower * ones(size(artifactsTrial(:,1))),'filled');
% scatter((artifactsTrial(:,2)-1)/Fs,threshLower * ones(size(artifactsTrial(:,2))),'filled');
% plot([t(1) t(end)],[threshUpper threshUpper]);
% plot([t(1) t(end)],[threshLower threshLower]);

end