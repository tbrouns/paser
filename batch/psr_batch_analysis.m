function psr_batch_analysis(subjectName,loadPathRoot,savePathRoot,expTypes,cfg)

if (nargin < 1); subjectName  = []; end
if (nargin < 2); loadPathRoot = []; end
if (nargin < 3); savePathRoot = []; end
if (nargin < 4); expTypes     = []; end

folderNames = dir([loadPathRoot 'R*']);
folderNames = char(folderNames.name);

%% sort files

nFolders = size(folderNames,1);
folders  = cell(nFolders,1);
for iFolder = nFolders:-1:1
    foldername = lower(folderNames(iFolder,:));
    foldername = strtrim(foldername);
    k          = strfind(foldername,expTypes);
    if ~isempty(k) || isempty(expTypes)
        folders{iFolder,1} = folderNames(iFolder,:); 
    end
end

folders = folders(~cellfun('isempty',folders)); % remove empty cells
nFolders = length(folders);

for iFolder = 1:nFolders 
    foldername = strtrim(folders{iFolder});
    
    loadPath = [loadPathRoot, foldername];
    savePath = [savePathRoot, foldername];
    
    if (loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end %#ok
    if (savePath(end) ~= '\'); savePath = [savePath, '\']; end %#ok
    
    MATfiles  = dir([loadPath '\PSR_' subjectName '*.mat']);
    nFiles    = size(MATfiles,1);
    filenames = char(MATfiles.name);
    if (isempty(filenames)); continue; end
    
    if (cfg.analysis.run)
        psr_wrapper_analysis(filenames,loadPath,savePath,cfg.analysis);
    end
    
    for iFile = 1:nFiles
        close all
        filename = filenames(iFile,:);
        fpath = [loadPath filename];
        load(fpath);
        [~,filename,~] = fileparts(filename);
                
        if (cfg.cluster.run)
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
        
        if (cfg.manual.run)
            if (~isfield(spikes.clusters.metrics,'labels'))
                labels = psr_manual_labelling(spikes,metadata,parameters,freq);
                for iClust = 1:length(labels); spikes.clusters.metrics(iClust).labels = labels(iClust); end
                save(fpath,'spikes','-append');
            end
        end
    end
end

end