function psr_batch_analysis(cfg)

% PSR_BATCH_ANALYSIS - Batch data analysis function
% 
% Syntax:  psr_batch_analysis(cfg)
%
% Inputs:
%    cfg - See toolbox README on how to set these parameters.
%
% Outputs:
%    One or more MAT files. See toolbox README for further details.
%
% See also: PSR_WRAPPER_FUNCTION

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (isempty_field(cfg,'cfg.subject'));      cfg.subject      = [];    end % Which subject to load data from
if (isempty_field(cfg,'cfg.loadpath'));     cfg.loadpath     = [];    end % Where to load the processed data
if (isempty_field(cfg,'cfg.savepath'));     cfg.savepath     = [];    end % Where to save the analysed  data
if (isempty_field(cfg,'cfg.pattern'));      cfg.pattern      = [];    end % Pattern to look for in foldernames
if (isempty_field(cfg,'cfg.analysis.run')); cfg.analysis.run = false; end % Whether to run the analysis script
if (isempty_field(cfg,'cfg.plot.quality')); cfg.plot.quality = false; end % Plot spike sorting quality control figures
if (isempty_field(cfg,'cfg.plot.merges'));  cfg.plot.merges  = false; end % Plot the cluster merges

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
    
    % Grab the processed data files
    
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
        psr_wrapper_function(cfg.analysis); % Wrapper function for calling custom analysis function
    end
    
    % Quality control: visualize unit quality and merges

    if (cfg.plot.quality || cfg.plot.merges)
        
        for iFile = 1:nFiles
            
            % Initialize
            
            close all
            spikes     = [];
            metadata   = [];
            parameters = [];
            
            % Load data 
            
            filename = filenames(iFile,:);
            filepath = [loadPath filename];
            load(filepath);
            [~,filename,~] = fileparts(filename);

            savePathFigs = [savePath 'clusters\'];
            
            if (~isempty_field(spikes,'spikes.spiketimes'))

                savePathQuality = [savePathFigs 'quality\'];
                savePathMerges  = [savePathFigs  'merges\'];
                
                if (cfg.plot.quality); [~,~,~] = mkdir(savePathQuality); psr_sst_plot_multiple(spikes,metadata,parameters,savePathQuality,filename); end % Plot the unit quality         
                if (cfg.plot.merges);  [~,~,~] = mkdir(savePathMerges);  psr_sst_plot_merges  (spikes,         parameters,savePathMerges, filename); end % Plot the merges 
            end
        end
    end
end

end