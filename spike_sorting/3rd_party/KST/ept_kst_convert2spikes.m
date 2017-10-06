function spikes = ept_kst_convert2spikes(rez,data,parameters)

Fs = parameters.Fs;
window_samples = round(Fs * parameters.spikes.window_size / 1000);
samples_before = round(Fs * parameters.spikes.cross_time  / 1000);
samples_after  = round(Fs * parameters.spikes.max_jitter  / 1000) + window_samples - (1 + samples_before);
win = -samples_before:samples_after;

% fid = fopen(rez.ops.fbinary, 'r');
% NchanTOT = rez.ops.NchanTOT;
% data = fread(fid, [NchanTOT inf], '*int16');
% fclose(fid);

nsamples = size(data,2);

% extract info from rez
spikeTimes     = rez.st3(:,1);
spikeTemplates = rez.st3(:,2);
if (size(rez.st3,2) >= 5)
    spikeClusters = 1+rez.st3(:,5);
else
    spikeClusters = spikeTemplates;
end

% get raw data around spiketimes
% waves = NaN(size(dat,1),numel(win),numel(spikeTimes));

times = bsxfun(@plus,spikeTimes,win);
times(times < 1) = 1;
times(times > nsamples) = nsamples;
times = times';
times = times(:);
waves = data(:,times);
waves = permute(waves,[3 2 1]);
waves = reshape(waves,length(win),[],parameters.general.nelectrodes);
waves = permute(waves,[2 1 3]);

% for i = 1:length(spikeTimes)
%     spkwin = spikeTimes(i) + win;
%     waves(:,:,i) = dat(:,spkwin);
% end

% find channel index with maximum amplitude template for each cluster
% [Currently not saved to 'spikes']
peakChannel = zeros(size(spikeClusters));
uClusters = unique(spikeClusters);
for c = 1:length(uClusters)
    clust     = uClusters(c);
    I         = spikeClusters == clust;
    templates =  unique(spikeTemplates(I));
    
    t = squeeze(range((rez.dWU(:,:,templates)),1));
    m = max(max(t));
    if any(size(t) ==1)
        chidx = find(t == m);
    else
        [chidx, ~] = find(t == m);
    end
    peakChannel(I) = chidx;
    
end

data                        = single(data');
spikes                      = [];
spikes.assigns              = single(spikeClusters');
spikes.spiketimes           = single(spikeTimes' / Fs);
spikes.waveforms            = single(waves);
spikes.info.detect.dur      = nsamples / Fs;
spikes.info.detect.thresh   = -parameters.spikes.thresh * (median(abs(data)) / 0.6745);
spikes.info.detect.stds     = std(data);
spikes.info.detect.cov      = get_covs({data}, window_samples);