function psr_sst_plot_multiple(spikes,metadata,parameters,freq,savePath,filename)

% PSR_SST_PLOT_MULTIPLE - Visualization of spike cluster obtained by PASER.
% This function creates a multi-plot of various visualization methods to
% access spike sorting quality and saves it as a PNG image.
%
% Syntax:  psr_sst_plot_multiple(spikes,metadata,parameters,freq,savePath,filename)
%
% Inputs:
%    spikes     - See README
%    metadata   - See README
%    parameters - See README
%    freq       - See README
%    savePath   - Directory to save output images in
%    filename   - Base filename of each image
% 
% Outputs:
%    PNG images in output folder specified by "savePath".
%
% See also: PSR_BATCH_VISUALIZATION

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

if (nargin < 4); freq     = []; end
if (nargin < 5); savePath = []; end % Save in current working directory
if (nargin < 6); filename = []; end

if (~isempty(filename)); [~,filename,~] = fileparts(filename); 
else,                    filename = 'Spikes'; 
end

nClust = size(spikes.clusters.vars,2);

% Convert and set necessary parameters

spikes = psr_sst_display_parameters(spikes);
spikes.unwrapped_times          = spikes.spiketimes;
spikes.params.detect.ref_period = parameters.spikes.ref_period; 
spikes.params.detect.shadow     = 0.5 * parameters.spikes.window_size;
spikes.info.kmeans.assigns      = spikes.assigns; 
numclusts                       = max(spikes.info.kmeans.assigns); % perhaps redundant? see nClust
cmap                            = jetm(numclusts);
spikes.info.kmeans.colors       = cmap(randperm(numclusts),:);

% Plot

fig = figure; set(gcf,'position',get(0,'screensize'));

for iClust = 1 : nClust
    
    clusterID = spikes.clusters.vars(iClust).id;
        
    if length(find(spikes.assigns == clusterID)) > parameters.cluster.min_spikes 
        figure(fig);
        clf
        spikes = psr_sst_filter_zscore(spikes,parameters,'delete');
        
        subaxis(2,4,[1 2],'Margin',0.05,'Padding',0); psr_sst_plot_waveforms  (spikes,clusterID,parameters);
        subaxis(2,4,    3,'Margin',0.05,'Padding',0); plot_detection_criterion(spikes,clusterID);
        subaxis(2,4,    4,'Margin',0.05,'Padding',0); psr_sst_plot_xcorr      (spikes,clusterID);
        subaxis(2,4,[5 6],'Margin',0.05,'Padding',0); psr_sst_plot_stability  (spikes,clusterID,freq,metadata,parameters);
        subaxis(2,4,    7,'Margin',0.05,'Padding',0); psr_sst_plot_poisson    (spikes,clusterID,parameters);
        subaxis(2,4,    8,'Margin',0.05,'Padding',0); plot_isi                (spikes,clusterID);
        
        export_fig([savePath filename '_C'  num2str(clusterID,'%02d')]);
    end
end

end

%------------- END OF CODE --------------