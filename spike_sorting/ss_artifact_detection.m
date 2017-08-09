function spikes = ss_artifact_detection(data,spikes,index_start)

window_size = (spikes.params.artifact_p2ptime / 1000) * spikes.params.Fs;  % in samples
num_chans   = size(data,2);
Fs          = spikes.params.Fs;

for ichan = 1:num_chans
    
    data_chan = data(:,ichan)';
        
    % Find peaks and troughs
    [pks_min,loc_min] = findpeaks(double(-data_chan));
    [pks_max,loc_max] = findpeaks(double( data_chan));
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
    
    d = d + 1;
    loc_all = [];
    
    for i = 1:d
        j   = N - (d - i);
        dn  = abs(loc_1(i:j) - loc_2);
        p2p = abs(pks_1(i:j) - pks_2);
        loc = mean([loc_1(i:j);loc_2],1); % between-peak locations
        id  = (dn <= window_size);
        p2p = p2p(id);
        loc = loc(id);
        stdev   = median(abs(data_chan)) / 0.6745; % median absolute deviation
        acutoff = spikes.params.artifact_thresh * stdev;
        loc_all = [loc_all; round(loc(p2p > acutoff))']; %#ok, artifact locations
    end
    
    spiketimes = (loc_all + index_start) / Fs;
    
    % Detect artifacts
    
    if (isfield(spikes,'artifacts')); spikes.artifacts = [spikes.artifacts; spiketimes]; % append
    else                              spikes.artifacts = spiketimes;
    end   
    
end

spikes.artifacts = single(spikes.artifacts)';

end