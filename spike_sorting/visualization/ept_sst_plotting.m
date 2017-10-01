function ept_sst_plotting(spikes,clusters,parameters,freq,which,savePath,filename)

if (nargin < 4); freq = []; end;
if (nargin < 5); which = 1:size(clusters.vars,2); end
if (nargin < 6); savePath = []; end % save in current working directory
if (nargin < 7); filename = []; end

if (~isempty(filename)); [~,filename,~] = fileparts(filename); 
else                     filename = 'Spikes'; 
end

clusters.vars = clusters.vars(which);

% Only plot single units % DO THRESHOLDING SOMEWHERE ELSE

% if (strcmp(type,'single'))
%     tf = ~strcmp({clusters.vars(:).unit},{'single'});
%     % Only plot significant amplitude units
%     tf = ~cell2mat({clusters.vars.flag});
%     clusters.vars(tf) = [];
% end

nclus = size(clusters.vars,2);

% Convert and set necessary parameters

spikes = ept_sst_display_parameters(spikes);
spikes.unwrapped_times          = spikes.spiketimes;
spikes.params.detect.ref_period = parameters.spikes.ref_period; 
spikes.info.kmeans.assigns      = spikes.assigns; 
numclusts                       = max(spikes.info.kmeans.assigns); % redundant? see nclus
cmap                            = jetm(numclusts);
spikes.info.kmeans.colors       = cmap(randperm(numclusts),:);

% Plot

fig = figure; set(gcf,'position',get(0,'screensize'));

for i = 1 : nclus
    
    icluster = clusters.vars(i).id;
    
    if length(find(spikes.assigns == icluster)) > parameters.cluster.size_min 
        figure(fig);
        clf
        subplot(2,3,1); plot_waveforms          (spikes,icluster);
        subplot(2,3,2); plot_residuals          (spikes,icluster);
        subplot(2,3,3); plot_distances          (spikes,icluster);
        subplot(2,3,4); plot_detection_criterion(spikes,icluster);
        subplot(2,3,5); plot_isi                (spikes,icluster);
        subplot(2,3,6); plot_stability          (spikes,icluster,freq,parameters);
        export_fig([savePath filename '_C'  num2str(icluster,'%03d')]);
    end
end

end