function batch_visualization(subject,loadPath_root,savePath_root)

if (nargin < 1); loadPath_root = []; end
if (nargin < 2); savePath_root = []; end

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
    else           folders_passive{iFolder,1} = folder_names(iFolder,:);
    end
end

folders_active  = folders_active (~cellfun('isempty',folders_active )); % remove empty cells
folders_passive = folders_passive(~cellfun('isempty',folders_passive)); % remove empty cells

numfolders_active  = length(folders_active);
numfolders_passive = length(folders_passive);

for iFolder = 1:numfolders_passive
    foldername   = folders_passive{iFolder};
        
    loadPath = [loadPath_root, foldername];
    savePath = [savePath_root, foldername];
    
    if (loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end
    if (savePath(end) ~= '\'); savePath = [savePath, '\']; end
    
    folder_names  = dir(loadPath);
    folder_names  = char(folder_names.name);
    numfolders    = length(folder_names(:,1));
    baseline      = [];
    
    MATfiles  = dir([loadPath '\Spikes_' subject '*.mat']);
    numfiles  = size(MATfiles,1);
    filenames = char(MATfiles.name);
    for iFile = 1:numfiles
        filename = filenames(iFile,:);
        load([loadPath filename]);
        spikes.params.cluster_method = 'fmm';
        [~,filename,~] = fileparts(filename);
        ss_visualization(spikes,clusters,'all',savePath,filename)
        if (metadata.stimulus == 0); baseline = freq;
        else                         ept_timefrequency_plotting(freq,parameters,savePath,filename,baseline);
        end
    end
    batch_analysis(MATfiles,loadPath,savePath);
end
end