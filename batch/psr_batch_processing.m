function psr_batch_processing(parameters)

% PSR_BATCH_PROCESSING - Batch data processing function for PASER.
%
% Syntax:  psr_batch_processing(parameters)
%
% Inputs:
%    parameters - See toolbox README on how to set these parameters.
%
% Outputs:
%    One or more MAT files. See toolbox README for further details.
%
% Subfunctions: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

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
if (~isfield(parameters,'txtfile'));    parameters.txtfile    = [];           end

folderNames = dir(parameters.loadPath); % Locate all files and folders in directory
folderNames = folderNames(~ismember({folderNames.name},{'.','..'}));
folderNames = folderNames([folderNames.isdir]); % Only keep folders
folderNames = char(folderNames.name); % Convert to character array

if (isempty(folderNames)); disp('No folders in directory'); return; end

OVERWRITE  = true;

if (~isempty(parameters.txtfile)) % Read sessions from text file
    filename = parameters.txtfile;
    fileID = fopen(filename);
    foldersTemp = textscan(fileID,'%s %s');
    nSessions = size(foldersTemp,2);
    for iSession = 1:nSessions
        folders(:,iSession) = foldersTemp{iSession};
    end
    fclose(fileID);
else
    %% Sort files
    
    nFolders = length(folderNames(:,1));
    
    nTypes = length(parameters.patterns);
    folderTypes = cell(nTypes+1,1);
    for iFolder = 1:nFolders
        folderName = lower(folderNames(iFolder,:));
        folderName = strtrim(folderName);
        k = zeros(nTypes,1);
        for iType = 1:nTypes
            k(iType) = ~isempty(strfind(folderName,lower(parameters.patterns{iType})));
        end
        k = find(k,1) + 1; % Find first match
        if (isempty(k)); k = 1; end % No pattern found
        folderTypes{k}{iFolder,1} = folderNames(iFolder,:);
    end
    
    nTypes = length(folderTypes);
    for iType = 1:nTypes % Remove empty cells
        type = folderTypes{iType};
        if (iscell(type))
            folderTypes{iType} = type(~cellfun('isempty',type));
        end
    end
    
    % Extract folders that match given type
    
    if (strcmp(parameters.type,'all') || isempty(parameters.patterns))
        folders = [];
        for iType = 1:nTypes
            folders = [folders;folderTypes{iType}];
        end
    else
        type = parameters.type;
        k = strfind(type,parameters.patterns) + 1;
        folders = folderTypes{k};
    end
end

% Determine which folders to process

switch parameters.process
    case 'given'
        nFolders = length(parameters.folders);
        if (nFolders > 0)
            foldersTemp = cell(nFolders,1);
            for iFolder = 1:nFolders
                k = strfind(folders,parameters.folders{iFolder});
                k = find(~cellfun(@isempty,k),1);
                if (~isempty(k))
                    foldersTemp{iFolder} = folders{k};
                end
            end
            folders = foldersTemp;
        else
           disp('The "process" field has been set to "given", but no folders were chosen in "folders" field.');
           return;
        end
    case 'new'
        OVERWRITE = false;
end

% Process data

if (isempty(folders)); disp('No folders found.'); return; end

nFolders = size(folders,1);
for iFolder = 1:nFolders
    
    nSessions = size(folders,2);
    
    parameters.session      = []; % session names
    parameters.sessionIndex = []; % note session for each trial
    parameters.loadPathSub  = []; % where each trial is loaded from
    
    FILES_FOUND = false;
    
    for iSession = 1:nSessions % sessions to process together
        
        folderName = strtrim(folders{iFolder,iSession});
        
        % Locate relevant files in subfolders
        loadPathSub = dir([parameters.loadPath, folderName '\**\*.' parameters.extension]);
        loadPathSub =    char(loadPathSub.folder);
        loadPathSub =  unique(loadPathSub,'rows');
        if (~isempty(loadPathSub)); FILES_FOUND = true; end
        loadPathSub = cellstr(loadPathSub);
        
        parameters.session{iSession} = folderName;
        parameters.loadPathSub  = [parameters.loadPathSub;loadPathSub];
        parameters.sessionIndex = [parameters.sessionIndex;iSession*ones(size(loadPathSub))];
        
        if (iSession == 1); parameters.savePathSub = folderName;
        else,               parameters.savePathSub = [parameters.savePathSub '-' folderName];
        end
    end
    
    if (FILES_FOUND)
        parameters.savePathSub = [parameters.savePath, parameters.savePathSub];
        
        if (~isempty(parameters.loadPathSub))
            [~,~,~] = mkdir(parameters.savePathSub);
            
            if (~OVERWRITE) % Check if 'spikes' files exist in save directory
                filesSpikes = dir([parameters.savePathSub '\PSR_*.mat']);
                filesTemp   = dir([parameters.savePathSub   '\Temp_*.mat']);
                filesSpikes = char(filesSpikes.name);
                filesTemp   = char(filesTemp.name);
                if size(filesSpikes,1) > 0 && size(filesTemp,1) == 0; continue; end
            end
            
            psr_wrapper(parameters);
        end
    end
end

%------------- END OF CODE --------------