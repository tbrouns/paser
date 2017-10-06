function batch_processing(subject,loadPath_root,savePath_root,expType) 

if (nargin < 1); subject = []; end
if (nargin < 2); loadPath_root = []; end
if (nargin < 3); savePath_root = []; end
if (nargin < 4); expType = 'all'; end

ptrn = 'passive';

folder_names = dir([loadPath_root 'R*']);
folder_names = char(folder_names.name);

if (isempty(folder_names)); disp('No folders in directory'); return; end

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

if     (strcmp(expType,'all'));     folders = [folders_active;folders_passive];
elseif (strcmp(expType,'active'));  folders = folders_active;
elseif (strcmp(expType,'passive')); folders = folders_passive;
end

numfolders = length(folders);

for iFolder = 7:numfolders
    foldername1   = strtrim(folders{iFolder});
    loadPath_sub1 = [loadPath_root, foldername1];
    savePath_sub1 = [savePath_root, foldername1];
    folder_names  = dir(loadPath_sub1);
    folder_names  = char(folder_names.name);
    numfolders    = length(folder_names(:,1));
    
    loadPath_sub2 = [];
    itr           = 1;
    
    for jFolder = 1:numfolders % across all trials 
        foldername2 = strtrim(folder_names(jFolder,:));
        if (~isempty(strfind(foldername2,'_')) || ~isempty(strfind(foldername2,'EPhys'))) % either different trials for passive or EPhys folder for active
            path = [loadPath_sub1 '\' foldername2]; 
            if (isdir(path)) % make sure we are storing folders and not files
                loadPath_sub2{itr} = path; %#ok
                itr = itr + 1;
            end
        end
    end
    
    if (~isempty(loadPath_sub2))
        ept_wrapper(subject,foldername1,loadPath_sub2,savePath_sub1);
    end
end