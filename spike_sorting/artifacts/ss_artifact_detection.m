function spikes = ss_artifact_detection(data,spikes,index_start)

method = 1;
visualize_0 = 0;
visualize_1 = 0;
visualize_2 = 0;
close all;

window_size = (spikes.params.mfa_p2ptime / 1000) * spikes.params.Fs;  % in samples
Fs          =  spikes.params.Fs;
tStart      = (index_start - 1) / Fs;
nlength     = size(data,1);
nchan       = size(data,2);

data_mean = mean(bsxfun(@rdivide,bsxfun(@minus,data,mean(data)),std(data)),2);
data_mean = data_mean';

if (visualize_1)
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(1:length(data_mean),data_mean);
    pause(1);
end

% Find peaks and troughs
[pks_min,loc_min] = findpeaks(double(-data_mean));
[pks_max,loc_max] = findpeaks(double( data_mean));
pks_min           = -pks_min;

% Calculate p2p amplitudes of adjacent peak+trough pair

Nmin = size(pks_min,2);
Nmax = size(pks_max,2);

if (Nmin > Nmax)
    pks_min = pks_min(1:Nmax);
    loc_min = loc_min(1:Nmax);
elseif (Nmax > Nmin)
    pks_max = pks_max(1:Nmin);
    loc_max = loc_max(1:Nmin);
end

if (loc_min(1) < loc_max(1));
    loc_1 = loc_min;
    loc_2 = loc_max;
    pks_1 = pks_min;
    pks_2 = pks_max;
else
    loc_1 = loc_max;
    loc_2 = loc_min;
    pks_1 = pks_max;
    pks_2 = pks_min;
end

% Detect putative artifacts based on amplitudes

loc_all = [];
p2p_all = [];

stdev   = median(abs(data_mean)) / 0.6745; % median absolute deviation
acutoff = spikes.params.mfa_thresh * stdev;

for i = 0:1
    
    dn  =   abs(loc_1(1+i:end)-loc_2(1:end-i));
    loc = mean([loc_1(1+i:end);loc_2(1:end-i)],1); % between-peak locations
    %     loc = loc_1(i:j); % align on one of the peaks
    
    p1  = pks_1(1+i:end);
    p2  = pks_2(1:end-i);
    
    if (sum(p1) < sum(p2)); pmin = p1; pmax = p2;
    else                    pmin = p2; pmax = p1;
    end
    
    % Filter based on p2p time
    id   = (dn <= window_size);
    pmax = pmax(id);
    pmin = pmin(id);
    loc  =  loc(id);
    
    p2p = pmax - pmin;
    
    % Filter based on max and min amplitudes
    id      = pmax > acutoff & pmin < -acutoff;
    loc_all = [loc_all, loc(id)]; %#ok, artifact locations
    p2p_all = [p2p_all, p2p(id)]; %#ok
end

spikeTimes = loc_all / Fs;

criterion = 0;
[artifacts,artifactTimesAll,~,nsamples] = ss_artifact_correlation(spikes,data_mean,spikeTimes,criterion);

data      = data';
waveforms = data(:,artifacts);
waveforms = permute(waveforms,[2 3 1]);
waveforms = reshape(waveforms,nsamples,[],nchan);
waveforms = single(permute(waveforms,[2 1 3]));
waveforms = reshape(waveforms,[],size(waveforms,2) * size(waveforms,3));

% Artifact sorting

assigns = ss_dictionary_learning(spikes,waveforms,spikes.params.mfa_fmm_p);

clustIDs = unique(assigns);
nclusts  = length(clustIDs);

MFA_half_width = round(0.5 * nsamples);

if (visualize_0);
    artifact_visualization(waveforms,assigns,nclusts);
    artifact_visualization_2(data_mean,artifactTimesAll,assigns,nclusts,Fs,MFA_half_width)
end

clear waveforms

% Find periodicity in clusters

% Fr = round(Fs / nsamples);

period    =  0;
valMax    =  0;
offset    =  0;
signalMax = [];
indexMax  =  0;
artifactTimesClustersOld  = cell(nclusts,1);

for iclust = 1:nclusts
    
    spikeID = (assigns == clustIDs(iclust));
    artifactTimes = artifactTimesAll(spikeID);
    
    n = round(Fs * artifactTimes)';
    artifactTimesClustersOld{iclust} = n;
    signal = createSignal(n,MFA_half_width,nlength);
    
    %     signal = resample(signal,Fr,Fs);
    %     signal(signal > 0) = 1;
    %     signal(signal < 0) = 0;
    
    % [pxx,f]  = periodogram(signal,[],[],Fs);
    % plot(f,pxx); % plot frequency power spectrum
    % [pxx,f]  = periodogram(     r,[],[],Fr);
    
    % Find periodicities in signal
    [r,lags] = xcorr(signal,'coeff'); % autocorrelation
    maxlag = round(Fs * (1 / spikes.params.mfa_freq_min));
    id = abs(lags) < maxlag;
    lNew = lags(id);
    rNew = r(id);
    stdev = median(abs(rNew)) / 0.6745; % median absolute deviation
    acutoff = spikes.params.mfa_std_corr * stdev;
    
    if (method)
        
        % METHOD I
        
        id1 = lNew >  Fs * (1 / spikes.params.mfa_freq_max);
        id2 = lNew < -Fs * (1 / spikes.params.mfa_freq_max);
        rNew = rNew(id1) + flipud(rNew(id2));
        lNew = lNew(id1);
        M = movsum(rNew,nsamples);
        [val,indexMaximum] = max(M);
        
        if (visualize_2)
            figure
            bar(lNew,rNew); % plot cross-correlations
            ylim([0 1.0]);
            
            %             figure
            %             bar(lNew,M);
        end
        
        if (val > valMax);
            
            M = M / nsamples;
            M = M < acutoff;
            I = find(M);
            I = I - indexMaximum;
            I1 = min(I(I >= 0));
            I2 = max(I(I <= 0));
            
            offset = round(0.5 * (I1 - I2));
            period = lNew(round(indexMaximum + mean([I1;I2]))) / Fs; % stimulus period (seconds)
            valMax = val;
            indexMax = iclust;
        end
        
    else
        
        % METHOD II
        
        if (acutoff < spikes.params.mfa_thr_corr); acutoff = spikes.params.mfa_thr_corr; end
        rThresh = rNew; % thresholded signal
        rThresh(rThresh < acutoff) = 0;
        
        [~,locs] = findpeaks(rThresh);
        locs = lNew(locs);
        locs = locs(locs > Fs * (1 / spikes.params.mfa_freq_max));
        
        % Find strongest periodicity
        
        for i = 1:length(locs)
            v = (locs(i):locs(i):lNew(end));
            width = floor(0.5 * locs(i) * spikes.params.mfa_freq_off);
            v = bsxfun(@plus,v,(-width:width)');
            v = v(:);
            v = v(v > 0);
            v = v(v <= lNew(end));
            v = [-flipud(v) ; v]; %#ok
            v = v + 1 - lNew(1);
            val = sum(rThresh(v)); % test
            if (val > valMax);
                period = locs(i) / Fs; % stimulus period (seconds)
                signalMax = signal;
                valMax = val;
            end
        end
    end
end

if (period < (1 / spikes.params.mfa_freq_max) || period > (1 / spikes.params.mfa_freq_min));
    return; % no strong periodicity
end

if (offset < MFA_half_width); offset = MFA_half_width; end

%% NEW METHOD

if (method)
    
    artifactTimesNew = [];
    
    for iclust = 1:nclusts
        
        artifactTimes = artifactTimesClustersOld{iclust};
        signal        = createSignal(artifactTimes,offset,nlength);
        nperiod       = round(Fs * period);

        N = zeros(length(artifactTimes),1);

        for n = 1:spikes.params.mfa_range

            t1 = artifactTimes + n * nperiod;
            t2 = artifactTimes - n * nperiod;

            t1 = bsxfun(@plus,t1,-offset:offset);
            t2 = bsxfun(@plus,t2,-offset:offset);

            t1(t1 <= 0) = 1;
            t1(t1 > nlength) = nlength;

            t2(t2 <= 0) = 1;
            t2(t2 > nlength) = nlength;

            N1 = sum(signal(t1),2);
            N2 = sum(signal(t2),2);
            N = N + (N1 > 0) + (N2 > 0);

        end

        % Adjust threshold for artifacts near start/end boundaries
        range = spikes.params.mfa_range;
        W1 = floor(         artifactTimes  / nperiod); 
        W2 = floor((nlength-artifactTimes) / nperiod);
        W  = min([W1 W2],[],2);
        id = W >= range;
        W( id) = 2 * range;
        W(~id) = range + W(~id); 
        id = N > W * spikes.params.mfa_frac;

        if (sum(id) / length(id)) >= spikes.params.mfa_clus_min
            artifactTimes = artifactTimes(id);
            artifactTimesNew = [artifactTimesNew;artifactTimes]; %#ok
%             figure;
%             plot_spikes((artifactTimes/Fs)',iclust,nclusts);
        end
    end
    
    if (~isempty(artifactTimesNew))
        MFA_times = sort((artifactTimesNew/Fs)+tStart);
        spikes.artifacts        = single(MFA_times);
        spikes.artifacts_period = period;
        spikes.artifacts_waves  = getWaveforms(data_mean,artifactTimesNew,MFA_half_width);
    end
    
else
    
    %% OLD METHOD
    
    % Create artificial signal to extract correct MFAs
    
    nstep = round(Fs*period);
    n     = (1:nstep:nlength)';
    signalBlock = createSignal(n,offset,nlength);
    
    % Find alignment point of artificial signal with regular signal
    % Acquire artifact times based on artificial signal
    
    [r,lags] = xcorr(signalMax,signalBlock,'coeff');
    [~,I]    = max(r);
    lagDiff  = lags(I);
    
    clear lags
    
    % Re-create artifical signal based on calculated lag
    
    if (lagDiff > 0); signalBlock = signalBlock( 1:end-lagDiff);
    else              signalBlock = signalBlock(-lagDiff+1:end);
    end
    
    signals = find(signalBlock);
    I = find(diff(signals) > 1);
    I = [0;I];
    I = round(I(2:end) - 0.5 * diff(I));
    signals      = signals(I);
    signalPrefix = signals(1)-nstep:-nstep:0;
    signalPrefix = (fliplr(signalPrefix))';
    signalSuffix = (signals(end)+nstep:nstep:nlength)';
    n = [signalPrefix;signals;signalSuffix];
    signalBlock = createSignal(n,offset,nlength);
    
    % Find MFAs based on all original artifact timings and created artificial
    % signal
    
    artifactTimesOld = [];
    for iclust = 1:nclusts
        artifactTimesOld = [artifactTimesOld;artifactTimesClustersOld{iclust}]; %#ok
    end
    
    signal = createSignal(artifactTimesOld,MFA_half_width,nlength);
    signalNew = signal .* signalBlock;
    
    artifactTimes = find(signalNew);
    I = find(diff(artifactTimes) > 1);
    if (~isempty(I));
        I = [0;I];
        I = round(I(2:end) - 0.5 * diff(I));
        artifactTimes = artifactTimes(I);
        artifactTimesNew = [];
        artifactTimesCal = [];
        nartifacts = length(artifactTimesOld);
        for iartifact = 1:nartifacts
            t = artifactTimesOld(iartifact);
            [M,I] = min(abs(t - artifactTimes));
            if (M <= offset)
                artifactTimesNew = [artifactTimesNew;t]; %#ok
                artifactTimesCal = [artifactTimesCal;artifactTimes(I)]; %#ok
            end
        end
        
        [artifactTimesNew,I] = sort(artifactTimesNew);
        artifactTimesCal = artifactTimesCal(I);
        
        while true % remove multiple artifacts at same position
            id = find(diff(artifactTimesNew) <= 2 * offset, 1);
            if (~isempty(id))
                t1 = artifactTimesNew(id);
                t2 = artifactTimesNew(id + 1);
                c1 = artifactTimesCal(id);
                c2 = artifactTimesCal(id + 1);
                [~,I] = max([abs(t1-c1);abs(t2-c2)]);
                artifactTimesNew(id+I-1) = [];
                artifactTimesCal(id+I-1) = [];
            else
                break;
            end
        end
        
        if (~isempty(artifactTimesNew))
            MFA_times = (artifactTimesNew/Fs)+tStart;
            spikes.artifacts        = single(sort(MFA_times));
            spikes.artifacts_period = period;
            spikes.artifacts_waves  = getWaveforms(data_mean,artifactTimesNew,MFA_half_width);
        end
    end
end
end

function signal = createSignal(n,nsamples,nlength)
n = bsxfun(@plus,n,-nsamples:nsamples);
n = n(:);
n = n(n > 0 & n <= nlength);
signal    = zeros(nlength,1);
signal(n) = 1;
end

function p2p = getP2Ps(data,times,width,Fs)

nlength = length(data);
n = round(Fs * times);
n = bsxfun(@plus,n,-width:width);
n(n <= 0) = 1;
n(n > nlength) = nlength;
data = data(n');
p2p  = max(data) - min(data);

end

function artifact_visualization(waveforms,assigns,nclusts)

for iclust = 1:nclusts
    
    figure
    
    id = assigns == iclust;
    waves = waveforms(id,:);
    
    cmap    = hot(64);
    [n,x,y] = histxt(waves);
    h       = imagesc(x,y,n);
    colormap(cmap);
    set(gca,'Color', cmap(1,:) );
    xlabel('Sample #');
    ylabel('Voltage');
    set(gca,'YDir','normal')
    
end

end

function artifact_visualization_2(dataAll,timesAll,assigns,nclusts,Fs,width)

nlength = length(dataAll);

for iclust = 1:nclusts
    
    id = (assigns == iclust);
    times = timesAll(id);
    times = round(Fs * times)';
    times = bsxfun(@plus,times,-width:width);
    times(times <= 0) = 1;
    times(times > nlength) = nlength;
    data = dataAll(times');
    figure;plot(1:size(data,1),data);
    
end

end

function data = getWaveforms(data,times,width)

nlength = length(data);
times = bsxfun(@plus,times,-width:width);
times(times <= 0) = 1;
times(times > nlength) = nlength;
data = data(times);

end

function artifact_visualization_3(dataAll,times,width)

nlength = length(dataAll);
times = bsxfun(@plus,times,-width:width);
times(times <= 0) = 1;
times(times > nlength) = nlength;
data = dataAll(times');
figure;plot(1:size(data,1),data);

end

function plot_spikes(spike_times,clustID,nclusts) 
    % Draw spikes (T is 1xN matrix with N spikes at times t)
    height = 0.1;
    plot([spike_times;spike_times],[ones(size(spike_times))*(clustID-height).';ones(size(spike_times))*(clustID+height).'],'k');
    ylim([0 nclusts]);
    xlabel('Time [s]');
end