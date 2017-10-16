function psr_batch_processing(parameters)

%% Check input

if (~isfield(parameters,'subject'));    parameters.subject    = [];           end
if (~isfield(parameters,'loadPath'));   parameters.loadPath   = [];           end
if (~isfield(parameters,'savePath'));   parameters.savePath   = [];           end
if (~isfield(parameters,'configPath')); parameters.configPath = [];           end
if (~isfield(parameters,'patterns'));   parameters.patterns   = [];           end
if (~isfield(parameters,'type'));       parameters.type       = 'all';        end
if (~isfield(parameters,'folders'));    parameters.folders    = [];           end
if (~isfield(parameters,'process'));    parameters.process    = 'all';        end
if (~isfield(parameters,'extension'));  parameters.extension  = 'continuous'; end

folder_names = dir(parameters.loadPath); % Locate all files and folders in directory
folder_names = folder_names(~ismember({folder_names.name},{'.','..'}));
folder_names = folder_names([folder_names.isdir]); % Only keep folders
folder_names = char(folder_names.name); % Convert to character array

if (isempty(folder_names)); disp('No folders in directory'); return; end

%% Sort files

numfolders = length(folder_names(:,1));

ntypes = length(parameters.patterns);
folder_types = cell(ntypes+1,1);
for iFolder = 1:numfolders
    foldername = lower(folder_names(iFolder,:));
    foldername = strtrim(foldername);
    k = zeros(ntypes,1);
    for iType = 1:ntypes
        k(iType) = ~isempty(strfind(foldername,lower(parameters.patterns{iType})));
    end
    k = find(k,1) + 1; % Find first match
    if (isempty(k)); k = 1; end % No pattern found
    folder_types{k}{iFolder,1} = folder_names(iFolder,:);
end

ntypes = length(folder_types);
for iType = 1:ntypes % Remove empty cells
    folder_types{iType} = folder_types{iType}(~cellfun('isempty',folder_types{iType}));
end

% Extract folders that match given type

if (strcmp(parameters.type,'all') || isempty(parameters.patterns))
    folders = [];
    for iType = 1:ntypes
        folders = [folders;folder_types{iType}];
    end
else
    type = parameters.type;
    k = strfind(type,parameters.patterns) + 1;
    folders = folder_types{k};
end

% Determine which folders to process

numfolders = length(parameters.folders);
OVERWRITE  = true;
switch parameters.process
    case 'given'
        foldersTemp = cell(numfolders,1);
        for iFolder = 1:numfolders
            k = strfind(folders,parameters.folders{iFolder});
            k = find(~cellfun(@isempty,k),1);
            if (~isempty(k))
                foldersTemp{iFolder} = folders{k};
            end
        end
        folders = foldersTemp;
    case 'new'
        OVERWRITE = false;
end

% Process data

if (isempty(folders)); disp('No folders found.'); return; end

numfolders = length(folders);
for iFolder = 1:numfolders
    foldername1   = strtrim(folders{iFolder});
    loadPath_sub1 = [parameters.loadPath, foldername1];
    savePath_sub1 = [parameters.savePath, foldername1];
    
    % Locate relevant files in subfolders
    loadPath_sub2 = dir([loadPath_sub1 '\**\*.' parameters.extension]);
    loadPath_sub2 = char(loadPath_sub2.folder); 
    loadPath_sub2 = unique(loadPath_sub2,'rows');
    loadPath_sub2 = cellstr(loadPath_sub2);
    
    parameters.session     = foldername1;
    parameters.loadPathSub = loadPath_sub2;
    parameters.savePathSub = savePath_sub1;
    
    if (~isempty(parameters.loadPathSub))
        [~,~,~] = mkdir(parameters.savePathSub);
        
        if (~OVERWRITE) % Check if 'spikes' files exist in save directory
            filesSpikes = dir([parameters.savePathSub '\Spikes_*.mat']);
            filesSpikes = char(filesSpikes.name);
            filesTemp   = dir([parameters.savePathSub   '\Temp_*.mat']);
            filesTemp   = char(filesTemp.name);
            if size(filesSpikes,1) > 0 && size(filesTemp,1) == 0; continue; end
        end
        
        psr_wrapper(parameters);
    end
end