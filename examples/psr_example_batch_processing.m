function psr_example_batch_processing(parameters)

% See "Quick Start" in README file on how to call this function

if (nargin < 1)
    psr_show_warning({...
        'ERROR: need to specify input paths'},...
        true);
    return;
end

% Set some optional flags 
PROCESS   = true; if (isfield(parameters,'processing')); PROCESS   = parameters.processing; end
VISUALIZE = true; if (isfield(parameters,'visualize'));  VISUALIZE = parameters.visualize;  end
ANALYSE   = true; if (isfield(parameters,'analyse'));    ANALYSE   = parameters.analyse;    end

% Set a couple variables
subjectID = 'SubjectID'; % Name of subject
savePath  = parameters.savePath; % Where we save the processed data

% Path of current file
fpath = mfilename('fullpath');
fpath = strsplit(fpath,'\');
fpath = cell2mat(join(fpath(1:end-1),'\'));
fpath = [fpath,'\'];

%% Data processing

if (PROCESS)
    
    stimPath = [fpath,'psr_example_stimulus']; % Define stimulus detection filepath
    
    parameters.stimPath     = stimPath;       % Function to detect stimulus onsets
    parameters.configPath   = [];             % Where the parameters are loaded from
    parameters.subject      = subjectID;      % Name of subject used in output MAT filename (no spaces)
    parameters.nelectrodes  = 4;              % Number of electrodes per polytrode (e.g. tetrode: 4)
    parameters.extension    = 'continuous';   % File extension of raw data files
    parameters.rawpattern   = 'CH';           % Pattern to look for in data files
    parameters.blockpattern = [];             % Used to differentiate between blocks within session
    parameters.stimpatterns = [];             % Which session type to process
    parameters.process      = 'all';          % Which specific sessions to process ('all', 'given' or 'from')
    parameters.folders      = [];             % Sessions that you wish to process, if 'given' or 'from' is chosen above (string cell array)
    parameters.overwrite    = false;          % Whether to overwrite data from existing processed sessions
    parameters.txtfile      = [];             % Folders to process given in text file
    
    psr_batch_processing(parameters); % Process raw data files
    
end

%% Visualize some of the processed data

if (VISUALIZE)
    
    % Units we want to do quality control for
    id_folder = 2;
    id_probe  = 1;
    
    % Which folder to load data from
    folders = dir([savePath '/' parameters.subject '*']);
    folders = char(folders.name);
    folder  = strtrim(folders(id_folder,:));
    
    % Which file to load
    files = dir([parameters.savePath '/' folder '/PSR*']);
    files = char(files.name);
    file  = strtrim(files(id_probe,:));
    
    % Load the processed data file
    load([parameters.savePath '/' folder '/' file]);
    
    % Figure file name 
    filename = [folder '_P' num2str(id_probe)];
    
    % Plot and save the quality control figures
    psr_sst_plot_multiple(spikes,metadata,parameters,savePath,filename);
    close all
    
end

%% Analyse the data

if (ANALYSE)
    
    cfg = [];
    cfg.loadpath = savePath;  % Location of output files from LFP and spike processing
    cfg.savepath = savePath;  % Where to save the output analysis files. We use the same location here.
    cfg.subject  = subjectID; % Subject name
    
    cfg.analysis.fpath      = [fpath,'psr_example_analysis']; % Location of the analysis function
    cfg.analysis.run        = true; % Whether we want to run the analysis
    cfg.analysis.configPath = [];   % Path to script containing custom analysis parameters (not used here)

    psr_batch_analysis(cfg);
    close all
    
end

end