function psr_batch_visualization(subjectName,loadPathRoot,savePathRoot,expType,cfg)

if (nargin < 1); subjectName  = [];    end
if (nargin < 2); loadPathRoot = [];    end
if (nargin < 3); savePathRoot = [];    end
if (nargin < 4); expType      = 'all'; end

ptrn = 'passive';

folderNames = dir([loadPathRoot 'R*']);
folderNames = char(folderNames.name);

%% sort files

nFolders = length(folderNames(:,1));

foldersPassive = cell(1,1);
foldersActive  = cell(1,1);
for iFolder = 1:nFolders
    foldername = lower(folderNames(iFolder,:));
    foldername = strtrim(foldername);
    k          = strfind(foldername,ptrn);
    if isempty(k); foldersActive {iFolder,1} = folderNames(iFolder,:);
    else,          foldersPassive{iFolder,1} = folderNames(iFolder,:);
    end
end

foldersActive  = foldersActive (~cellfun('isempty',foldersActive )); % remove empty cells
foldersPassive = foldersPassive(~cellfun('isempty',foldersPassive)); % remove empty cells

if     (strcmp(expType,'all'));     folders = [foldersActive;foldersPassive];
elseif (strcmp(expType,'active'));  folders = foldersActive;
elseif (strcmp(expType,'passive')); folders = foldersPassive;
end

nFolders = length(folders);

for iFolder = 11:nFolders
    foldername = strtrim(folders{iFolder});
    
    loadPath = [loadPathRoot, foldername];
    savePath = [savePathRoot, foldername];
    
    if (loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end %#ok
    if (savePath(end) ~= '\'); savePath = [savePath, '\']; end %#ok
    
    MATfiles  = dir([loadPath '\Spikes_' subjectName '*.mat']);
    nFiles    = size(MATfiles,1);
    filenames = char(MATfiles.name);
    if (isempty(filenames)); continue; end
    
    if (cfg.analysis.flag)
        savePathFigures = [savePath 'analysis\'];
        [~,~,~] = mkdir(savePathFigures);
        psr_batch_analysis(filenames,loadPath,loadPath,savePathFigures,cfg.analysis);
    end
    
    for iFile = 1:nFiles
        close all
        filename = filenames(iFile,:);
        fpath = [loadPath filename];
        load(fpath);
        [~,filename,~] = fileparts(filename);
        run(parameters.general.configPath); % TEMP
        if (cfg.cluster)
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
        end
        
        if (cfg.manual)
            if (~isfield(spikes.clusters.metrics,'labels'))
                labels = psr_manual_labelling(spikes,metadata,parameters,freq);
                for iClust = 1:length(labels); spikes.clusters.metrics(iClust).labels = labels(iClust); end
                save(fpath,'spikes','-append');
            end
        end
    end
end

end