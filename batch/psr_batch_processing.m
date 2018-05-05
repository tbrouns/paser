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
if (~isfield(parameters,'folders'));    parameters.folders    = [];           end
if (~isfield(parameters,'process'));    parameters.process    = 'all';        end
if (~isfield(parameters,'extension'));  parameters.extension  = 'continuous'; end
if (~isfield(parameters,'txtfile'));    parameters.txtfile    = [];           end
if (~isfield(parameters,'overwrite'));  parameters.overwrite  = false;        end

folderNames = dir(parameters.loadPath); % Locate all files and folders in directory
folderNames = folderNames(~ismember({folderNames.name},{'.','..'}));
folderNames = folderNames([folderNames.isdir]); % Only keep folders
folderNames = char(folderNames.name); % Convert to character array

if (isempty(folderNames)); disp('No folders in directory'); return; end

if (~isempty(parameters.txtfile)) % Read sessions from text file
    filename = parameters.txtfile;
    
    str = '%s';
    while true
        fileID = fopen(filename);
        foldersTemp = textscan(fileID,str);
        fclose(fileID);
        x = foldersTemp{end};
        I = cellfun(@isempty,x);
        if (all(I))
            foldersTemp = foldersTemp(1:end-1);
            break;
        else
            str = [str ' %s'];
        end
        
    end
    
    folders = cell(0,0);
    nSessions = size(foldersTemp,2);
    for iSession = 1:nSessions
        foldersNew = foldersTemp{iSession};
        nfolders = size(folders,1);
        if (length(foldersNew) < nfolders)
            foldersNew{nfolders} = [];
        end
        
        folders(:,iSession) = foldersNew;
    end
else
    %% Sort files
    
    patterns = [];
    nFolders = length(folderNames(:,1));
    folderTypes = cell(0);
    for iFolder = 1:nFolders
        folderNameRaw = folderNames(iFolder,:);
        folderNameRaw = strtrim(folderNameRaw);
        folderName = lower(folderNameRaw);
        folderName = split(folderName,'_');
        folderName = folderName(2:end);
        if (~isempty(folderName)); folderName = folderName{end};
        else,                      folderName = '';
        end
        k = strcmp(patterns,folderName);
        if (~k)
            patterns{end+1} = folderName;
            k = length(patterns);
            folderTypes{k,1} = [];
        end
        folderTypes{k,1}{end+1,1} = folderNameRaw;
    end
    
    % Extract folders that match given type
    
    if (isempty(parameters.patterns))
        folders = [];
        for iType = 1:nTypes
            folders = [folders;folderTypes{iType}];
        end
    else
        k = false(1,length(folderTypes));
        for iType = 1:length(parameters.patterns)
            k = k | strcmpi(patterns,parameters.patterns{iType});
        end
        folders = folderTypes(k);
    end
    
    folders = cat(1,folders{:});
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
                [k,~] = ind2sub(folders,k);
                if (~isempty(k))
                    foldersTemp(iFolder,:) = folders(k,:);
                end
            end
            folders = foldersTemp;
        else
            disp('The "process" field has been set to "given", but no folders were chosen in "folders" field.');
            return;
        end
    case 'from'
        k = strfind(folders,parameters.folders{1});
        k = find(~cellfun(@isempty,k),1);
        [k,~] = ind2sub(size(folders),k);
        folders = folders(k:end,:);
end

% Process data

if (isempty(folders)); disp('No folders found.'); return; end

nFolders = size(folders,1);
for iFolder = 1:nFolders
    
    nSessions = size(folders,2);
    
    parameters.session      = []; % session names
    parameters.sessionIndex = []; % note session for each trial
    parameters.loadPathSub  = []; % where each trial is loaded from
    parameters.savePathSub  = [];
    
    for iSession = 1:nSessions % sessions to process together
        folderName = folders{iFolder,iSession};
        if (~isempty(folderName))
            folderName = strtrim(folderName);
            % Locate relevant files in subfolders
            loadPathSub = dir([parameters.loadPath, folderName '\**\*.' parameters.extension]);
            loadPathSub =    char(loadPathSub.folder);
            loadPathSub =  unique(loadPathSub,'rows');
            if (~isempty(loadPathSub))
                loadPathSub = cellstr(loadPathSub);
                parameters.session{iSession} = folderName;
                parameters.loadPathSub  = [parameters.loadPathSub;loadPathSub];
                parameters.sessionIndex = [parameters.sessionIndex;iSession*ones(size(loadPathSub))];
                if (isempty(parameters.savePathSub)); parameters.savePathSub = [parameters.subject     '-' folderName];
                else,                                 parameters.savePathSub = [parameters.savePathSub '-' folderName];
                end
            end
        end
    end
    
    if (~isempty(parameters.loadPathSub))
        parameters.savePathSub = [parameters.savePath, parameters.savePathSub];
        [~,~,~] = mkdir(parameters.savePathSub);
        if (~parameters.overwrite) % Check if 'spikes' files exist in save directory
            filesSpikes = dir([parameters.savePathSub  '\PSR_*.mat']);
            filesTemp   = dir([parameters.savePathSub '\Temp_*.mat']);
            filesSpikes = char(filesSpikes.name);
            filesTemp   = char(filesTemp.name);
            if size(filesSpikes,1) > 0 && size(filesTemp,1) == 0; continue; end
        end
        
        psr_wrapper(parameters);
    end
end

%------------- END OF CODE --------------