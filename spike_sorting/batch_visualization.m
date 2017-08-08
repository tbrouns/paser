function batch_visualization(loadPath_root,savePath_root) 

if (nargin < 1); loadPath_root = []; end
if (nargin < 2); savePath_root = []; end

% loadPath_root = 'D:\Data\electrophysiology\YZ05 (Mouse 34)\';
% savePath_root = 'D:\Data\electrophysiology\YZ05 (Mouse 34) - Figures';

ptrn = 'passive';

folder_names = dir('R*');
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
    foldername1   = folders_passive{iFolder};
    loadPath_sub1 = [loadPath_root, foldername1];
    savePath_sub1 = [savePath_root, foldername1];
    folder_names  = dir(loadPath_sub1);
    folder_names  = char(folder_names.name);
    numfolders    = length(folder_names(:,1));
    
    MATfiles  = dir([loadPath_sub1 '\Spikes_YZ*.mat']);
    numfiles  = size(MATfiles,1);
    filenames = char(MATfiles.name);
    for iFile = 1:numfiles
       load([loadPath_sub1 '\' filenames(iFile,:)]);
       spikes.params.cluster_method = 'fmm';
       ss_visualization(spikes,clusters,'all',savePath_sub1,filenames(iFile,:))
    end
end
end