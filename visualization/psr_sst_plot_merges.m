function psr_sst_plot_merges(spikes,parameters,savePath,filename)

if (nargin < 3); savePath = []; end % save in current working directory
if (nargin < 4); filename = []; end

if (~isempty(filename)); [~,filename,~] = fileparts(filename); 
else,                    filename = 'Spikes'; 
end

plotNum = 6;
plotIds = fliplr([1,2,3,4,5,8]);

parameters = psr_load_parameters(parameters,'display');
parameters.display.metrics = false;

assigns          = spikes.assigns;
assignsPrior     = spikes.assigns_prior;
spikes.assigns   = spikes.assigns_prior;
spikes.waveforms = psr_int16_to_single(spikes.waveforms,parameters);

clustPriorId = unique(assignsPrior);
clustMergeId = zeros(size(clustPriorId),class(clustPriorId));

nClustsPrior = length(clustPriorId);
for iClust = 1:nClustsPrior
    clustMergeId(iClust) = assigns(find(assignsPrior == clustPriorId(iClust),1));
end

clustIds = unique(clustMergeId);
nClusts  = length(clustIds);
for iClust = 1:nClusts
    clustMain = clustIds(iClust);
    clustSub = clustPriorId(clustMergeId == clustMain);
    clustSub = clustSub(clustSub ~= clustMain);
    
    nClustSub = length(clustSub);
    
    if (nClustSub > 0)
        
        nFigs = ceil(nClustSub / plotNum);
        for iFig = 1:nFigs
            figure; set(gcf,'position',get(0,'screensize'));
            
            iStart = plotNum * (iFig - 1) + 1;
            iEnd   = plotNum * iFig;
            
            if (iEnd > length(clustSub)); clustFig = clustSub(iStart:end);
            else,                         clustFig = clustSub(iStart:iEnd);
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
            
            subaxis(2,4,[6 7]);
            psr_sst_plot_waveforms(spikes,clustMain,parameters);
            suptitle(['Merges of Cluster #' num2str(clustMain) ' (' num2str(iFig) ' of ' num2str(nFigs) ')']);
            
            if (nFigs > 1); str = ['_' num2str(iFig)];
            else,           str = [];
            end
            
            export_fig([savePath filename '_C'  num2str(clustMain,'%02d') str]);
        end
    end
end

end