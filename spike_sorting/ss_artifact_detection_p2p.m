function spikes = ss_artifact_detection(data,spikes,index_start)

window_size = (spikes.params.artifact_p2ptime / 1000) * spikes.params.Fs;  % in samples

% Calculate z-score
data_mean = mean(data,1);
zscore    = data_mean - mean(data_mean) / std(data_mean);
clear data_mean;

% Find peaks and troughs
[pks_min,loc_min] = findpeaks(double(-zscore));
[pks_max,loc_max] = findpeaks(double( zscore));
pks_min           = -pks_min;

% Calculate p2p amplitudes of adjacent peak+trough pair

Nmin = size(pks_min,2);
Nmax = size(pks_max,2);

if (Nmin >= Nmax)
    pks_1 = pks_min;
    pks_2 = pks_max;
    loc_1 = loc_min;
    loc_2 = loc_max;
    d     = Nmin - Nmax;
    N     = Nmin;
else
    pks_1 = pks_max;
    pks_2 = pks_min;
    loc_1 = loc_max;
    loc_2 = loc_min;
    d     = Nmax - Nmin;
    N     = Nmax;
end

d   = d + 1;
loc_all = [];

for i = 1:d
    j   = N - (d - i);
    dn  = abs(loc_1(i:j) - loc_2);
    p2p = abs(pks_1(i:j) - pks_2);
    loc = mean([loc_1(i:j);loc_2],1); % between-peak locations
    id  = (dn <= window_size);
    p2p = p2p(id);
    loc = loc(id);
    stdev   = median(abs(zscore)) / 0.6745; % median absolute deviation
    acutoff = spikes.params.artifact_thresh * stdev;
    loc_all = [loc_all; round(loc(p2p > acutoff))']; % artifact locations
end

% Detect artifacts

if (isfield(spikes,'artifacts')); spikes.artifacts = [spikes.artifacts; loc_all + index_start];
else                              spikes.artifacts = loc_all + index_start;
end

spikes.artifacts = single(spikes.artifacts);

end