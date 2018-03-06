function data = psr_ms_denoise_raw(data,parameters,Fs)

% Deal with large signal transitions to reduce ringing artifacts after filtering

sWinSlope    = round(Fs * parameters.ms.denoise.raw.win.slope    / 1000);
sWinArtifact = round(Fs * parameters.ms.denoise.raw.win.artifact / 1000);
sWinPadding  = round(Fs * parameters.ms.denoise.raw.win.padding  / 1000);
sWinPulse    = round(Fs * parameters.ms.denoise.raw.win.pulse    / 1000);

% First detect

derivative = data(1+sWinSlope:end) - data(1:end-sWinSlope);
sLength    = length(derivative);
thresh     = parameters.ms.denoise.raw.thresh * psr_mad(derivative);
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
artifacts(artifacts(:,2) > sLength, 2) = sLength;

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