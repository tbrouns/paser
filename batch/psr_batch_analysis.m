function psr_batch_analysis(cfg)

if (psr_isempty_field(cfg,'cfg.subject'));  cfg.subject  = []; end
if (psr_isempty_field(cfg,'cfg.loadpath')); cfg.loadpath = []; end
if (psr_isempty_field(cfg,'cfg.savepath')); cfg.savepath = []; end
if (psr_isempty_field(cfg,'cfg.pattern'));  cfg.pattern  = []; end

folderNames = dir([cfg.loadpath cfg.subject '*']);
folderNames = char(folderNames.name);

%% Sort files

nFolders = size(folderNames,1);
folders  = cell(nFolders,1);
for iFolder = nFolders:-1:1
    foldername = lower(folderNames(iFolder,:));
    foldername = strtrim(foldername);
    k          = strfind(foldername,cfg.pattern);
    if ~isempty(k) || isempty(cfg.pattern)
        folders{iFolder,1} = strtrim(folderNames(iFolder,:)); 
    end
end

folders = folders(~cellfun('isempty',folders)); % remove empty cells
nFolders = length(folders);

for iFolder = 1:nFolders
    foldername = folders{iFolder};
    
    loadPath = [cfg.loadpath, foldername];
    savePath = [cfg.savepath, foldername];
    
    if (loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end %#ok
    if (savePath(end) ~= '\'); savePath = [savePath, '\']; end %#ok
    
    MATfiles  = dir([loadPath '\PSR_' cfg.subject '*.mat']);
    nFiles    = size(MATfiles,1);
    filenames = char(MATfiles.name);
    if (isempty(filenames)); continue; end
    
    % LFP and spike analysis
    if (cfg.analysis.run)
        cfg.analysis.files    = filenames;
        cfg.analysis.loadpath = loadPath;
        cfg.analysis.savepath = savePath;
        psr_wrapper_function(cfg.analysis);
    end
        
    if (cfg.cluster.run) || (cfg.manual.run)
        
        for iFile = 1:nFiles
            close all
            filename = filenames(iFile,:);
            filepath = [loadPath filename];
            load(filepath);
            [~,filename,~] = fileparts(filename);

            % Visualize spike clusters
            if (cfg.cluster.run)
                savePathClusters = [savePath 'clusters\'];
                if (isfield(spikes,'spiketimes'))

                    savePathQuality = [savePathClusters 'quality\'];
                    [~,~,~] = mkdir(savePathQuality);
                    psr_sst_plot_multiple(spikes,metadata,parameters,savePathQuality,filename);

                    savePathMerges = [savePathClusters 'merges\'];
                    [~,~,~] = mkdir(savePathMerges);
                    psr_sst_plot_merges(spikes,parameters,savePathMerges,filename);
                end
            end

            % Spike cluster labelling (for development)
            if (cfg.manual.run)
                if (~isfield(spikes.clusters.metrics,'labels'))
                    labels = psr_manual_labelling(spikes,metadata,parameters,freq);
                    for iClust = 1:length(labels); spikes.clusters.metrics(iClust).labels = labels(iClust); end
                    save(filepath,'spikes','-append');
                end
            end
        end
    end
end

end