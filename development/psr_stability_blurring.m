function signal = psr_stability_blurring(spikes,signal,parameters)

clusterIDs = unique(spikes.assigns);
nClust     = length(clusterIDs);
nChan      = size(signal,1);
sLength    = size(signal,2);
precision  = 10^parameters.general.precision;
Fs         = parameters.Fs;

window_samples = round(Fs * parameters.spikes.window_size / 1000);
samples_hwidth = round(0.5 * window_samples);
win = (-samples_hwidth:samples_hwidth)';

for iClust = 1:nClust
    
    % Spike indices
    clusterID = clusterIDs(iClust);
    which     = spikes.assigns == clusterID;
    t = spikes.spiketimes(which);
    t = round(t * Fs + 1);
    t(t <= samples_hwidth)           = [];
    t(t >  sLength - samples_hwidth) = [];
    nspikes = length(t);
    t = bsxfun(@plus,t,win);
    t = t(:);
    
    % Extract waveforms
    w = psr_single(signal(:,t),parameters);
    w = permute(w,[3 2 1]);
    w = reshape(w,[],nspikes,nChan);
    w = permute(w,[2 1 3]);
    
    % Blurring
    wmean = mean(w,1); % Mean waveform of cluster
    wp = w(randperm(nspikes),:,:); % permutated waveforms
    wp = bsxfun(@minus,wp,wmean); % relative to mean waveform
    w = w + parameters.develop.gamma * wp; % blurred
    
    % Insert waveforms
    w = permute(w,[2 1 3]);
    w = reshape(w,1,[],nChan);
    w = permute(w,[3 2 1]); 
    
    signal(:,t) = int16(precision * w);
end

end