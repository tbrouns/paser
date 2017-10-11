function ept_sst_plotting(spikes,parameters,freq,which,savePath,filename,label)

if (nargin < 3); freq = []; end
if (nargin < 4); which = 1:size(spikes.clusters.vars,2); end
if (nargin < 5); savePath = []; end % save in current working directory
if (nargin < 6); filename = []; end
if (nargin < 7); label = []; end

if (~isempty(filename)); [~,filename,~] = fileparts(filename); 
else,                    filename = 'Spikes'; 
end

spikes.clusters.vars = spikes.clusters.vars(which);

% Only plot single units % DO THRESHOLDING SOMEWHERE ELSE

% if (strcmp(type,'single'))
%     tf = ~strcmp({spikes.clusters.vars(:).unit},{'single'});
%     % Only plot significant amplitude units
%     tf = ~cell2mat({spikes.clusters.vars.flag});
%     spikes.clusters.vars(tf) = [];
% end

nclus = size(spikes.clusters.vars,2);

% Convert and set necessary parameters

spikes = ept_sst_display_parameters(spikes);
spikes.unwrapped_times          = spikes.spiketimes;
spikes.params.detect.ref_period = parameters.spikes.ref_period; 
spikes.params.detect.shadow     = 0.5 * parameters.spikes.window_size;
spikes.info.kmeans.assigns      = spikes.assigns; 
numclusts                       = max(spikes.info.kmeans.assigns); % redundant? see nclus
cmap                            = jetm(numclusts);
spikes.info.kmeans.colors       = cmap(randperm(numclusts),:);

% Plot

fig = figure; set(gcf,'position',get(0,'screensize'));

for i = 1 : nclus
    
    icluster = spikes.clusters.vars(i).id;
    
    if length(find(spikes.assigns == icluster)) > parameters.cluster.min_spikes 
        figure(fig);
        clf
        spikes = ept_sst_filter_zscore(spikes,parameters,'delete');
        subplot(2,3,1); plot_waveforms          (spikes,icluster);
        subplot(2,3,2); plot_residuals          (spikes,icluster);
        subplot(2,3,3); plot_distances          (spikes,icluster);
        subplot(2,3,4); plot_detection_criterion(spikes,icluster);
        subplot(2,3,5); plot_isi                (spikes,icluster);
        subplot(2,3,6); plot_stability          (spikes,icluster,freq,parameters);
        export_fig([savePath filename '_C'  num2str(icluster,'%03d') label]);
    end
    
    % [x1,x2,~] = plot_fld(spikes, clusterID, clusterIDs(jcluster), 0); %
    % CURRENTLY NOT USED
end

end