function [artifactTimes_1,artifactTimes_2] = ss_mfa_combine(files,artifactTimesControl_1,artifactTimesControl_2)

% Grab artifacts across tetrodes

ntets    = length(files);
dur      = zeros(ntets,1);
off      = zeros(ntets,1);
periods  = zeros(ntets,1);
MFAtimesAll = cell(ntets,1);

iTetrodeMax = 0;
frac_max = 0;
for iTetrode = 1:ntets
    load(files{iTetrode});
    if isfield(spikes,'artifacts')
        periods(iTetrode) = spikes.artifacts_period;
        MFAtimesAll{iTetrode} = sort(spikes.artifacts);
        nmax = (spikes.info.detect.dur / spikes.artifacts_period);
        frac = length(spikes.artifacts) / nmax;
        if (frac > frac_max)
            iTetrodeMax = iTetrode;
            frac_max = frac;
        end
    end
    dur(iTetrode) = spikes.info.detect.dur;
    off(iTetrode) = spikes.params.mfa_freq_off;
end

dur = mode(dur);
off = mode(off);

artifactTimes = [];

if (iTetrodeMax ~= 0)
    
    period_main  = periods(iTetrodeMax);
    
    % create artificial signal
    
    MFAtimes = MFAtimesAll{iTetrodeMax}';
    MFAtimesPrefix = MFAtimes(1)-period_main:-period_main:0;
    MFAtimesPrefix = fliplr(MFAtimesPrefix);
    MFAtimesSuffix = MFAtimes(end)+period_main:period_main:dur;
    Ts = [MFAtimesPrefix MFAtimes MFAtimesSuffix]';
    
    while true
        dT = diff(Ts);
        dT = round(dT / period_main);
        id = find(dT > 1, 1);
        if (isempty(id)); break; end
        d = dT(id);
        Tm = linspace(Ts(id),Ts(id+1),d+1)';
        Ts = [Ts(1:id);Tm(2:end-1);Ts(id+1:end)];
    end
    
    % Grab artifact timings of clusters with similar main period
    
    MFAtimesNew = [];
    
    for iTetrode = 1:ntets
        if (periods(iTetrode) < period_main * (1 + off) && periods(iTetrode) > period_main * (1 - off))
            MFAtimesNew = [MFAtimesNew;MFAtimesAll{iTetrode}];
        end
    end
    
    MFAtimesNew = sort(MFAtimesNew);

    nartifacts = length(Ts);
        
    for iArtifact = 1:nartifacts
        T = Ts(iArtifact);
        times = MFAtimesNew(abs(MFAtimesNew - T) < period_main * off);
        if (~isempty(times));
            artifactTimes = [artifactTimes;mean(times)]; %#ok
        end
    end
end

% Couple with control artifacts 

offset_1 = findOffset(artifactTimes,artifactTimesControl_1);
offset_2 = findOffset(artifactTimes,artifactTimesControl_2);

if (~isnan(offset_1)); artifactTimes_1 = artifactTimesControl_1 + offset_1;
else                   artifactTimes_1 = artifactTimesControl_1;
end

if (~isnan(offset_2)); artifactTimes_2 = artifactTimesControl_2 + offset_2;
else                   artifactTimes_2 = artifactTimesControl_2;
end

end

function offset = findOffset(T,Tc)

[X,Y] = meshgrid(T,Tc);
[D,I] = min(abs(X-Y));
distances = [];
D2 = X - Y;
D_signed = zeros(1,length(I));

for i = 1:length(I)
   D_signed(i) = D2(I(i),i);
end

for i = 1:max(I)
    I1 = find(I == i);
    [~,I2] = min(D(I1));
    I2 = I1(I2);
    D2 = D_signed(I2);
    if (~isempty(D2))
        distances = [distances;D2]; %#ok
    end
end

dsorted = sort(abs(distances));
d = median(dsorted(1:round(0.5 * length(dsorted))));
I = abs(distances) < 2 * d;
offset = mean(distances(I));

end
