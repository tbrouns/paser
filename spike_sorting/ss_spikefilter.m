function [spikes,clusters] = ss_spikefilter(spikes,clusters)

clustID = cell2mat({clusters.vars.id});
flags   = cell2mat({clusters.vars.flag});
clustID = clustID(flags);

nclusts = length(clustID);

% Remove outliers based on absolute amplitude

nspikes  = size(spikes.waveforms,1);
nsamples = size(spikes.waveforms,2);
nchans   = size(spikes.waveforms,3);

waves  = reshape(spikes.waveforms,nspikes,nsamples*nchans);
id     = max(sign(spikes.params.outlier_abs) * waves,[],2) >= sign(spikes.params.outlier_abs) * spikes.params.outlier_abs;
spikes = ss_spike_removal(spikes,~id);

% Remove outliers based on z-score

for iclust = 1:nclusts
    [x,y,z,~] = plot_distances(spikes, clustID(iclust), 1, 0);
    y         = y / sum(y);
    [~,I1]    = max(y);
    I2        = find(y(I1:end) < spikes.params.outlier_chi,1);
    I         = I2 + I1 - 1;
    z_thresh  = x(I);
    if (~isempty(z_thresh))
        id        = z > z_thresh;
        which     = find(spikes.assigns == clustID(iclust));
        id        = which(id);
        spikes    = ss_spike_removal(spikes,id);
    end
end

% Remove outliers based on PC distance

for iclust = 1:nclusts
    
    which = find(spikes.assigns == clustID(iclust));
    x = spikes.waveforms(which,:) * spikes.info.pca.v(:,1);
    y = spikes.waveforms(which,:) * spikes.info.pca.v(:,2);
    z = spikes.waveforms(which,:) * spikes.info.pca.v(:,3);

    M  = [x,y,z];
    D  = mahal(M,M);
    sd = std(D);
    id = D > spikes.params.outlier_std * sd;
    id = which(id);
    
    spikes = ss_spike_removal(spikes,id);
    
end

% Deal with RPVs

spikes = ss_rpv_filter(spikes);

end

function spikes = ss_rpv_filter(spikes)

clusters = ss_clusterfeatures(spikes);

clustID = cell2mat({clusters.vars.id});
flags   = cell2mat({clusters.vars.flag});
clustID = clustID(flags);

nclusts = length(clustID);

for iclust = 1:nclusts

show = get_spike_indices(spikes, clustID(iclust));

x = spikes.waveforms(show,:) * spikes.info.pca.v(:,1);
y = spikes.waveforms(show,:) * spikes.info.pca.v(:,2);
z = spikes.waveforms(show,:) * spikes.info.pca.v(:,3);

PC_M    = [x,y,z];

spiketimes   = spikes.unwrapped_times(show);
rpvs         = diff(spiketimes) <= 0.001 * spikes.params.refractory_period;
rpvs         = [0,rpvs]; %#ok
rpvs         = find(rpvs);
id           = zeros(size(spiketimes));
id(rpvs)     = 1;
id(rpvs - 1) = 1;
id           = find(id);

num_rpvs = length(id);
itr = 1;
del = [];
n   = 0;

while (itr < num_rpvs)
    if (spiketimes(id(itr+1)) - spiketimes(id(itr)) <= 0.001 * spikes.params.refractory_period)
        v     = [id(itr);id(itr+1)];
        itrV  = [itr;itr+1];
        [~,I] = max(mahal(PC_M(v,:),PC_M));
        I1    = itrV(I);
        I2    = v(I);
        spiketimes(I2) = [];
        PC_M(I2,:)     = [];
        id(I1:end)     = id(I1:end) - 1;
        id(I1)         = [];
        del = [del;I2+n]; %#ok
        num_rpvs = num_rpvs - 1;
        itr = itr - 1;
        n   = n + 1;
    end
    
    itr = itr + 1;
end

id = show(del);
spikes.rpvs = length(id);
spikes = ss_spike_removal(spikes,id);
 
end

end