function psr_sst_plot_merges(spikes,parameters,savePath,filename)

if (nargin < 3); savePath = []; end % save in current working directory
if (nargin < 4); filename = []; end

if (~isempty(filename)); [~,filename,~] = fileparts(filename); 
else,                    filename = 'Spikes'; 
end

fig = figure; set(gcf,'position',get(0,'screensize'));

% Do conversion to avoid conversion for each separate plot
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

% Plot locations
plotNum = 6;
plotIds = fliplr([1,2,3,4,5,8]);

% Load display parameters
parameters = psr_parameters_load(parameters,'display');
parameters.display.metrics = false;

% Set pre-merge assigns and post-merge assigns
assigns          = spikes.assigns; % Post-merge assigns
assignsPrior     = spikes.assigns_prior; % Pre-merge assigns

clustPriorIDs = unique(assignsPrior); % All pre-merge cluster IDs
clustMergeIDs = zeros(size(clustPriorIDs),class(clustPriorIDs)); % Corresponding post-merge cluster IDs
 
for iClust = 1:length(clustPriorIDs)
    clustMergeIDs(iClust) = assigns(find(assignsPrior == clustPriorIDs(iClust),1));
end

clustIDs = unique(clustMergeIDs); % All post-merge cluster IDs
nClusts  = length(clustIDs);

for iClust = 1:nClusts
    
    clustMainID = clustIDs(iClust); % The main cluster that the sub-clusters were merged into
    clustSubIDs = clustPriorIDs(clustMergeIDs == clustMainID); % Clusters that formed new main cluster
    nClustSub = length(clustSubIDs);
    
    if (nClustSub > 0)
        
        nFigs = ceil(nClustSub / plotNum);
        for iFig = 1:nFigs
            
            figure(fig); clf;
            
            %% Plot pre-merge sub-clusters
            
            spikes.assigns = assignsPrior; % Use old assigns
            
            iStart = plotNum * (iFig - 1) + 1;
            iEnd   = plotNum * iFig;
            
            if (iEnd > length(clustSubIDs)); clustFig = clustSubIDs(iStart:end);
            else,                            clustFig = clustSubIDs(iStart:iEnd);
            end
            
            nPlots = length(clustFig);
            for iPlot = 1:nPlots
                
                subaxis(2,4,plotIds(iPlot),'PaddingTop',0.015);
                psr_sst_plot_waveforms(spikes,clustFig(iPlot),parameters);
                
                if (nPlots == 1) % Padding
                    h = subaxis(2,4,plotIds(iPlot + 1));
                    plot(0,0);
                    set(h,'Visible','off')
                end
            end
            
            %% Plot post-merge main cluster         
            
            spikes.assigns = assigns; % Use new assigns
            
            subaxis(2,4,[6 7]);
            psr_sst_plot_waveforms(spikes,clustMainID,parameters);
            suplabel(['Merges of Cluster #' num2str(clustMainID) ' (' num2str(iFig) ' of ' num2str(nFigs) ')']);
            
            if (nFigs > 1); str = ['_' num2str(iFig)];
            else,           str = [];
            end
            
            export_fig([savePath filename '_C'  num2str(clustMainID,'%02d') str]);
        end
    end
end

end