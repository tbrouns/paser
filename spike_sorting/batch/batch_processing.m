function batch_processing(subject,method,loadPath_root,savePath_root) 

if (nargin < 1); subject = []; end
if (nargin < 2 || isempty(method)); method = 'fmm'; end
if (nargin < 3); loadPath_root = []; end
if (nargin < 4); savePath_root = []; end

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

for iFolder = 2:numfolders_passive
    foldername1   = folders_passive{iFolder};
    loadPath_sub1 = [loadPath_root, foldername1];
    savePath_sub1 = [savePath_root, foldername1];
    folder_names  = dir(loadPath_sub1);
    folder_names  = char(folder_names.name);
    numfolders    = length(folder_names(:,1));
    
    loadPath_sub2 = [];
    itr           = 1;
    
    for jFolder = 1:numfolders % across all trials 
        foldername2 = strtrim(folder_names(jFolder,:));
        if ~isempty(strfind(foldername2,'_')) % check if underscore is present
            loadPath_sub2{itr} = [loadPath_sub1 '\' foldername2]; %#ok
            itr = itr + 1;
        end
    end
    
    if (~isempty(loadPath_sub2))
        ss_wrapper(subject,foldername1,loadPath_sub2,savePath_sub1,method);
    end
end