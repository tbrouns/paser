function [data,spikes] = ums_artifact_removal(data,spikes)

window_size     = (spikes.params.artifact_p2ptime / 1000) * spikes.params.Fs;  % in samples
artifact_length = floor(0.5 * (spikes.params.artifact_length / 1000) * spikes.params.Fs); % in samples

% Calculate z-score
data_mean = mean(data{1},2);
zscore    = data_mean - mean(data_mean) / std(data_mean);
clear data_mean;

% Find peaks and troughs
[pks_min,loc_min] = findpeaks(double(-zscore)); 
[pks_max,loc_max] = findpeaks(double( zscore));
pks_min           = -pks_min;

% Calculate p2p amplitudes of adjacent peak+trough pair

if (loc_max(1) < loc_min(1));    
    pks_max = pks_max(2:end);   
    loc_max = loc_max(2:end);
end

Nmin = size(pks_min,1);
Nmax = size(pks_max,1);

if (Nmin > Nmax); 
    pks_min = pks_min(1:Nmax); 
    loc_min = loc_min(1:Nmax);
elseif (Nmax > Nmin) % probably redundant
    pks_max = pks_max(1:Nmin); 
    loc_max = loc_max(1:Nmin);
end

dn  = loc_max - loc_min;
p2p = pks_max - pks_min;
loc = mean([loc_max,loc_min],2); % between-peak locations
id  = dn <= window_size;
p2p = p2p(id);
loc = loc(id);

% Detect artifacts

stdev           = median(abs(zscore)) / 0.6745; % median absolute deviation
acutoff         = spikes.params.artifact_thresh * stdev;
loc             = round(loc(p2p > acutoff)); % artifact locations
artifacts       = bsxfun(@plus,loc,-artifact_length:artifact_length);
artifacts       = artifacts(:);
artifacts       = unique(artifacts);
artifacts       = artifacts(artifacts > 0 & artifacts < size(zscore,1));
artifacts       = sort(artifacts);

data{1}(artifacts,:) = 0; % remove entirely
spikes.artifacts     = artifacts;

end