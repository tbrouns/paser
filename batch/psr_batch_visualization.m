function psr_batch_visualization(subject,loadPath_root,savePath_root,expType,CLUSTERS,ANALYSIS,PLOTTING)

if (nargin < 1); subject       = [];    end
if (nargin < 2); loadPath_root = [];    end
if (nargin < 3); savePath_root = [];    end
if (nargin < 4); expType       = 'all'; end
if (nargin < 5); CLUSTERS      = true;  end
if (nargin < 6); ANALYSIS      = true;  end
if (nargin < 7); PLOTTING      = true;  end % Plot figures during analysis

ptrn = 'passive';

folder_names = dir([loadPath_root 'R*']);
folder_names = char(folder_names.name);

%% sort files

numfolders = length(folder_names(:,1));

folders_passive = cell(1,1);
folders_active  = cell(1,1);
for iFolder = 1:numfolders
    foldername = lower(folder_names(iFolder,:));
    foldername = strtrim(foldername);
    k          = strfind(foldername,ptrn);
    if isempty(k); folders_active {iFolder,1} = folder_names(iFolder,:);
    else,          folders_passive{iFolder,1} = folder_names(iFolder,:);
    end
end

folders_active  = folders_active (~cellfun('isempty',folders_active )); % remove empty cells
folders_passive = folders_passive(~cellfun('isempty',folders_passive)); % remove empty cells

if     (strcmp(expType,'all'));     folders = [folders_active;folders_passive];
elseif (strcmp(expType,'active'));  folders = folders_active;
elseif (strcmp(expType,'passive')); folders = folders_passive;
end

numfolders = length(folders);

for iFolder = 1:numfolders
    foldername = strtrim(folders{iFolder});
    
    loadPath = [loadPath_root, foldername];
    savePath = [savePath_root, foldername];
    
    if (loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end %#ok
    if (savePath(end) ~= '\'); savePath = [savePath, '\']; end %#ok
            
    baseline  = [];
    MATfiles  = dir([loadPath '\Spikes_' subject '*.mat']);
    numfiles  = size(MATfiles,1);
    filenames = char(MATfiles.name);
    if (isempty(filenames)); continue; end
    
    if (ANALYSIS)
        savePathFigures = [savePath 'analysis\'];
        [~,~,~] = mkdir(savePathFigures);
        psr_batch_analysis(filenames,loadPath,loadPath,savePathFigures,PLOTTING);
    end
    
    if (CLUSTERS)
        for iFile = 1:numfiles
            close all
            filename = filenames(iFile,:);
            load([loadPath filename]);
            [~,filename,~] = fileparts(filename);
            savePathClusters = [savePath 'clusters\'];
            if (isfield(spikes,'spiketimes'))
                
                savePathQuality = [savePathClusters 'quality\'];
                [~,~,~] = mkdir(savePathQuality);
                psr_sst_plot_multiple(spikes,metadata,parameters,freq,savePathQuality,filename); 
                
                savePathMerges = [savePathClusters 'merges\'];
                [~,~,~] = mkdir(savePathMerges);
                psr_sst_plot_merges(spikes,parameters,savePathMerges,filename); 
                
                %             psr_sst_plot_clusters(spikes,clusters,savePath,filename);
            end
            
            %%% FIX LFP PLOTTING
%             if (~isempty(spikes.stimtimes{1}))
%                 if (metadata.stimulus == 0); psr_lfp_plotting(freq,parameters,savePath,filename); baseline = getBaseline(freq); % must be loaded first
%                 else,                        psr_lfp_plotting(freq,parameters,savePath,filename,baseline);
%                 end
%             end

        end
    end
end

end

function base = getBaseline(freq)
    base = freq.powspctrm;
    base = nanmean(base,3);
    base = repmat(base,1,1,size(freq.time,2));
end