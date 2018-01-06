function data = psr_artifact_diff(data,parameters,Fs)

% Deal with large signal transitions to reduce ringing artifacts after filtering

sWinSlope    = round(Fs * parameters.filter.diff.win_slope    / 1000);
sWinArtifact = round(Fs * parameters.filter.diff.win_artifact / 1000);
sWinPadding  = round(Fs * parameters.filter.diff.win_padding  / 1000);
sWinPulse    = round(Fs * parameters.filter.diff.win_pulse    / 1000);

% First detect

derivative = data(1+sWinSlope:end) - data(1:end-sWinSlope);
sLength    = length(derivative);
thresh     = parameters.filter.diff.thresh * psr_mad(derivative);
locThresh  = find(abs(derivative) > thresh);

% Group samples together to form artifacts

itr = 1;
N   = length(locThresh);
artifacts = [];

while itr < N
    locStart = locThresh(itr);
    width = 0;
    while width <= sWinArtifact
        locEnd = locThresh(itr+1);
        width = locEnd - locStart;
        itr = itr + 1;
        if (itr >= N); itr = itr + 1; break; end
    end
    artifacts = [artifacts;[locStart,locThresh(itr-1)]];
end

if isempty(artifacts); return; end

% Add padding

artifacts(:,1) = artifacts(:,1) - sWinPadding;
artifacts(:,2) = artifacts(:,2) + sWinSlope + sWinPadding;

artifacts(artifacts(:,1) < 1,       1) = 1;
artifacts(artifacts(:,2) > sLength, 1) = sLength;

% Find artifact pairs

iArtifact  = 1;
nArtifacts = size(artifacts,1);
artifactPairs = [];

while iArtifact < nArtifacts
    
    T1 = artifacts(iArtifact,:);
    id = find((artifacts(iArtifact+1:end,1) - T1(1)) <= sWinPulse);
    T2 = artifacts(iArtifact + id,:); % neighbouring artifacts
    
    % Find largest neighbouring artifact
    N       = size(T2,1);
    maxVals = zeros(N,1);
    for i = 1:N; maxVals(i) = max(abs(derivative(T2(i,1):T2(i,2)))); end
    [~,I] = max(maxVals);
    T2 = T2(I,:);
    
    if (~isempty(T2) && T2(1) > T1(2)) % Save
        artifactPairs = [artifactPairs;[T1,T2]];
        iArtifact = iArtifact + 1;
    end
    
    iArtifact = iArtifact + 1;
end

% Do signal subtraction
% 
% Calculate average signal
% 
% ndims = 3;
% nsamples = 20;
% nsections = 3;
% tRelax = 0.1; % [sec]
% 
% artifacts = artifactPairs(:,1);
% sWin = mean(diff(artifacts));
% sRlx = round(Fs * tRelax);
% d = sort(diff([artifactPairs(1:end-1,4),artifactPairs(2:end,1)],[],2));
% sEnd = round(0.9 * median(d(1:round(0.5 * length(d)))));
% 
% Segment data
% 
% del = artifacts + sWin > length(data);
% artifacts(del) = [];
% id = bsxfun(@plus,artifacts,0:sWin);
% dataSegmented = data(id); 
% 
% Initial dimensionality reduction through downsampling
% 
% N  = round(nsamples * [0.25,0.40,0.40]);
% S0 =  artifactPairs(:,1);
% S1 =  artifactPairs(:,2:3);
% S2 = [artifactPairs(:,4),artifactPairs(:,4)+sRlx];
% S3 = [S2(:,2)+1,S2(:,1)+sEnd];
% S1 = S1 - S0;
% S2 = S2 - S0;
% S3 = S3 - S0;
% 
% nArtifacts = length(artifacts);
% 
% Zcell    = cell(1,nsections); 
% Zcell{1} = NaN(nArtifacts,2*N(1)); % Initialize with buffer
% Zcell{2} = NaN(nArtifacts,2*N(2));
% Zcell{3} = NaN(nArtifacts,2*N(3));
% 
% for iArtifact = 1:nArtifacts
%     Y  = dataSegmented(iArtifact,:);
%     Y1 = Y(S1(iArtifact,1):S1(iArtifact,2));
%     Y2 = Y(S2(iArtifact,1):S2(iArtifact,2));
%     Y3 = Y(S3(iArtifact,1):S3(iArtifact,2));
%     f1 = round(Fs * N(1) / length(Y1));
%     f2 = round(Fs * N(2) / length(Y2));
%     f3 = round(Fs * N(3) / length(Y3));
%     y1 = resample(Y1,f1,Fs);
%     y2 = resample(Y2,f2,Fs);
%     y3 = resample(Y3,f3,Fs);
%     Zcell{1}(iArtifact,1:length(y1)) = y1;
%     Zcell{2}(iArtifact,1:length(y2)) = y2;
%     Zcell{3}(iArtifact,1:length(y3)) = y3;
% end
% 
% Z = [];
% for iSec = 1:nsections
%     id = find(sum(isnan(Zcell{iSec})),1);
%     if (id > 2)
%         Z = [Z,Zcell{iSec}(:,1:id-1)];
%     end
% end
% 
% Further dimensionality reduction through PCA
% 
% [~,PC] = pca(Z);
% PC = PC(:,1:ndims);
% mu = mean(PC);
% sd =  std(PC);
% secIDs = false(size(PC));
% 
% for iDim = 1:ndims
%     upperThresh = mu(iDim) + 2 * sd(iDim);
%     lowerThresh = mu(iDim) - 2 * sd(iDim);
%     secIDs(:,iDim) = PC(:,iDim) < upperThresh & PC(:,iDim) > lowerThresh;
% end
% 
% secIDs = (sum(secIDs,2) > ceil(0.5 * ndims));
% meanSection = median(dataSegmented(secIDs,:));
% 
% Downsample
% 
% artifactsNew = artifacts(secIDs,:);
% figure; hold on;
% plot(data);
% scatter(artifactsNew(:),zeros(size(artifactsNew(:))));
% 
% figure; 
% plot(meanSection);
% 
% figure;
% plot(mean(dataSegmented(secIDs,:)));
% 
% data = data - meanSection;

%% Adjust inter-artifact section

nArtifacts = size(artifactPairs,1);

for iArtifact = 1:nArtifacts
    
    T1 = artifactPairs(iArtifact,[1 2]);
    T2 = artifactPairs(iArtifact,[3 4]);
    
    % Extract between-peak window and adjust
    
    Y1 = data(T1(1));
    Y2 = data(T2(2));
    Ymean = 0.5 * (Y1 + Y2);
    
    I = T1(2):T2(1);
    pulseWindow = data(I);
    pulseWindow = pulseWindow - mean(pulseWindow) + Ymean;
    
    dY1 = Y1 - pulseWindow(1);
    dY2 = Y2 - pulseWindow(end);
    dY = linspace(dY1,dY2,length(pulseWindow))';
    pulseWindow = pulseWindow + dY;
    data(I) = pulseWindow;
    
    % Silence artifact
    
    data = silenceArtifact(data,T1(1):T1(2));
    data = silenceArtifact(data,T2(1):T2(2));
        
end

end

function data = silenceArtifact(data,artifactWindow)

% Interpolation

x  = [min(artifactWindow),max(artifactWindow)];
v  = data(x);
xq = artifactWindow;
vq = interp1(x,v,xq,'linear');
data(xq) = vq;

end