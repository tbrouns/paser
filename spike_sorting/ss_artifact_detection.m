function spikes = ss_artifact_detection(data,spikes,index_start)

window_size = (spikes.params.artifact_p2ptime / 1000) * spikes.params.Fs;  % in samples
Fs          = spikes.params.Fs;
tStart      = (index_start - 1) / Fs;
nchan       = size(data,2);

data_mean = mean(bsxfun(@rdivide,bsxfun(@minus,data,mean(data)),std(data)),2);
data_mean = data_mean';

offsetUpper = 1 / (spikes.params.artifact_freq - spikes.params.artifact_freqoff);
offsetLower = 1 / (spikes.params.artifact_freq + spikes.params.artifact_freqoff);

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

loc_all = [];
p2p_all = [];

stdev   = median(abs(data_mean)) / 0.6745; % median absolute deviation
acutoff = spikes.params.artifact_thresh * stdev;

% close all
% figure;
% hold on
% plot((1:size(data_mean,2)) / Fs, data_mean);
% plot((1:size(data_mean,2)) / Fs, acutoff*ones(size(data_mean)),'--r');
% plot((1:size(data_mean,2)) / Fs,-acutoff*ones(size(data_mean)),'--r');

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

% % Split data into sections based on stimulus frequency. Take maximum p2p
% % signal in each section.
%
% period    = 1 / (spikes.params.artifact_freq);
% nsections = ceil(2*(max(spikeTimes)/period));
% t1 = 0;
%
% spikeTimesNew = [];
%
% for isection = 1:nsections
%     t2 = t1 + period;
%     id = spikeTimes >= t1 & spikeTimes < t2;
%     times = spikeTimes(id);
%     p2p   = p2p_all(id);
%     if (~isempty(p2p))
%         [~,I] = max(p2p);
%         spikeTimesNew = [spikeTimesNew,times(I)]; %#ok
%     end
%     t1 = t1 + 0.5 * period;
% end
%
% spikeTimes = unique(spikeTimesNew);

criterion = 0;
[artifacts,artifactTimes,~,nsamples] = ss_artifact_correlation(spikes,data_mean,spikeTimes,criterion);

data = data';
waveforms = data(:,artifacts);
waveforms = permute(waveforms,[2 3 1]);
waveforms = reshape(waveforms,nsamples,[],nchan);
waveforms = single(permute(waveforms,[2 1 3]));
waveforms = reshape(waveforms,[],size(waveforms,2) * size(waveforms,3));

assigns = ss_dictionary_learning(spikes,waveforms,spikes.params.artifact_fmm_p);

nclusts  = length(unique(assigns));
clustIDs = unique(assigns);

MFA_times = [];

iclust = 1;
while (iclust <= nclusts)
    
    spikeID1 = (assigns == clustIDs(iclust));
    t1  = artifactTimes(spikeID1);
    dt = diff(t1);
    id = dt >= offsetLower & dt <= offsetUpper;
    f1  = sum(id) / length(id);
    
    jclust = iclust + 1;
    while (jclust <= nclusts)
        
        spikeID2 = (assigns == clustIDs(jclust));
        t2 = artifactTimes(spikeID2);
        dt = diff(t2);
        id = dt >= offsetLower & dt <= offsetUpper;
        f2 = sum(id) / length(id);
        
        T  = sort([t1,t2]);
        dt = diff(T);
        id = dt >= offsetLower & dt <= offsetUpper;
        F  = sum(id) / length(id);
        
        if (F > f1 && F > f2)
            f1 = F;
            t1 = T;
            assigns(spikeID2) = clustIDs(iclust);
            clustIDs(jclust)  = [];
            nclusts = nclusts - 1;
            jclust  = jclust  - 1;
        end
        
        jclust = jclust + 1;
    end
    
    if (f1 >= spikes.params.artifact_ratio && sum(spikeID1) > spikes.params.cluster_min)
        
        loc = find(assigns == clustIDs(iclust));
        
        % Remove waveforms with subthreshold period
        % Procedure is similar to RPV removal
        
        dt  = diff(t1);
        dt  = dt <= offsetLower;
        dt  = [0, dt];
        dt  = find(dt);
        
        id         = zeros(size(t1));
        id(dt)     = 1;
        id(dt - 1) = 1;
        id         = find(id);
        
        % Each "RPV" involves two or more spike. We remove enough spikes to resolve
        % the RPV, where we keep the spikes that have the largest p2p
        % amplitude 
        
        N = length(id);
        itr = 1;
        n   = 0;
        
        while (itr < N)
            if (t1(id(itr+1)) - t1(id(itr)) <= offsetLower)
                v     = [id(itr);id(itr+1)];
                itrV  = [itr;itr+1];
                
                % Find largest p2p amplitude (mean across channels)
                
                waves = waveforms(loc,:);
                waves = waves(v,:);
                amps  = zeros(length(v),1);
                for i = 1:length(v)
                   wave = waves(i,:);
                   wave = reshape(wave',nsamples,[]);
                   wave = mean(wave,2);
                   amps(i) = max(wave) - min(wave); 
                end
                
                [~,I] = min(amps); % remove smallest amplitude waveform                
                I1    = itrV(I);
                I2    = v(I);
                
                t1(I2)            = [];
                waveforms(I2,:)   = [];
                assigns(I2)       = [];
                artifactTimes(I2) = [];
                loc               = find(assigns == clustIDs(iclust));
                id(I1:end)        = id(I1:end) - 1;
                id(I1)            = [];
                
                N   =   N - 1;
                itr = itr - 1;
                n   =   n + 1;
            end
            
            itr = itr + 1;
        end
        
        MFA_times = [MFA_times,t1+tStart]; %#ok
        
%         figure;
%         subplot(1,2,1);
%         cmap    = hot(64);
%         [n,x,y] = histxt(single(waveforms(loc,:)));
%         h       = imagesc(x,y,n);
%         colormap(cmap);
%         set(gca,'Color', cmap(1,:) );
%         xlabel('Sample #');
%         ylabel('Voltage');
%         axis xy;
%         
%         subplot(1,2,2);
%         plot(1:size(waveforms,2),single(waveforms(loc,:)),'Color',[0.5 0.5 0.5])
%         xlabel('Sample #');
%         ylabel('Voltage');
    end
    
    iclust = iclust + 1;
    
end

spikes.artifacts = single(MFA_times);

disp(['Detected ' num2str(length(MFA_times)) ' magnetic field artifacts.']);

end
