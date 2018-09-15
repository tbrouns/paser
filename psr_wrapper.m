function psr_wrapper(parameters)

% PSR_WRAPPER - Wrapper for processing data with PASER.
% This function loads raw extracellular data, performs all data processing
% steps and saves the results to a MAT file. The processing pipeline
% consists of spike, local field potential and stimulus detection, as well
% as spike sorting.
%
% Syntax:  psr_wrapper(parameters)
%
% Inputs:
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    One or more MAT files. See toolbox README for further details.
%
% See also: PSR_BATCH_PROCESSING

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

%% Check input

subject  = []; if (isfield(parameters,'subject'));     subject  = parameters.subject;     end
session  = []; if (isfield(parameters,'session'));     session  = parameters.session;     end
loadPath = []; if (isfield(parameters,'loadPathSub')); loadPath = parameters.loadPathSub; end
savePath = []; if (isfield(parameters,'savePathSub')); savePath = parameters.savePathSub; end

%% Check if paths have correct format

nBlocks = length(loadPath);
for iBlock = 1:nBlocks
    if (loadPath{iBlock}(end) ~= '\' && ~isempty(loadPath{iBlock}))
        loadPath{iBlock} = [loadPath{iBlock}, '\']; %#ok
    end
end

if (savePath(end) ~= '\' && ~isempty(savePath)); savePath = [savePath, '\']; end

%% Convert initial parameters
temp.general = orderfields(parameters);
metadata     = temp.general;
parameters   = temp;
clear temp;

%% Move path field. Only relevant for example script
if (isfield(parameters.general,'path'))
    parameters.path    = parameters.general.path;
    parameters.general = rmfield(parameters.general,'path');
end

%% Load processing parameters
parameters = psr_parameters_load(parameters);

%% Constants
rawpattern  =    parameters.general.rawpattern;
precision   = 10^parameters.general.precision;
logfilePath = [savePath 'Temp_logfile.mat'];

%% Check if TEMP files exist in save folder

filesTemp = cell(0,0);
files  = dir([savePath 'Temp_Probe_*_Block_*.mat']);
files  = char(files.name);
nfiles = size(files,1);
for iFile = 1:nfiles
    filename = [savePath files(iFile,:)];
    iProbe = extractStringFromPath(filename,'Probe_(\d*)_');
    iBlock = extractStringFromPath(filename,'Block_(\d*).');
    filesTemp{iProbe,iBlock} = filename;
end
tempFilesFound = ~cellfun(@isempty,filesTemp);
logfileFound = psr_exist_in_file(logfilePath,'filesTemp');

%% Start data processing

if (~logfileFound || any(~tempFilesFound(:)) || isempty(tempFilesFound))
    
    %% Extract experimental conditions from folder name
    
    parameters.general.stimuli = cell(nBlocks,1);
    
    for iBlock = 1:nBlocks
        filepath = loadPath{iBlock};
        if (~isempty_field(parameters,'parameters.general.blockpattern'))
            parameters.general.stimuli{iBlock} = extractStringFromPath(filepath,['_(\d*)' parameters.general.blockpattern]);
        else
            string = regexp(filepath,'\','split');
            parameters.general.stimuli{iBlock} = string{end-1};
        end
    end
    
    %% Find raw data files
    
    filesUnsorted = cell(nBlocks,1);
    for iBlock = 1:nBlocks
        files = dir([loadPath{iBlock} '\*' rawpattern '*' parameters.general.extension]);
        if (size(files,1) == 0); return; end
        filesUnsorted{iBlock} = char(files.name);
    end
    
    %% Check raw data files and sort them in correct chronological order
    
    filesRawAll = cell(nBlocks,2);
    
    for iBlock = 1:nBlocks
        nFiles = length(filesUnsorted{iBlock}(:,1));
        filesRaw = cell(nFiles,2);
        for iFile = 1:nFiles
            filename   = filesUnsorted{iBlock}(iFile,:);
            filename   = strtrim(filename);
            itr        = strfind(filename,rawpattern) + length(rawpattern);
            [~,name,~] = fileparts(filename); % remove extension
            str        = name(itr:end); % take string after pattern
            k          = strfind(str,'_'); % check if underscore is present
            if (~isempty(k)); str = str(1:k-1); % take number between pattern and underscore
            else, k = length(str);
            end
            id = str2double(str); % convert to array index
            filesRaw{id,1} = filename;
            filesRaw{id,2} = [name(1:itr-1),name(itr+k:end)];
        end
        
        % Check consistency
        
        % Find filename that matches most other files
        matches = zeros(nFiles,1);
        for iFile = 1:nFiles
            filename = filesRaw{iFile,2}; % Compare file with others
            I = strcmp(filename,filesRaw(:,2));
            matches(iFile) = sum(I);
        end
        [~,I] = max(matches);
        I = strcmp(filesRaw{I,2},filesRaw(:,2));
        
        if (any(~I))
            str_1 = ['WARNING: raw data filenames are not consistent. Check data folder: "' loadPath{iBlock} '"'];
            if (parameters.process.warning)
                str   = cell2mat(join(cellstr(num2str(find(~I),' %d\n'))',','));
                str_2 = ['Ignoring channel(s): ' str '. Max. channel: ' num2str(find(I,1,'last'))];
            else
                str_2 = ['Processing terminated for session "' session '"'];
            end
            psr_show_warning({str_1,str_2},true);
            if (~parameters.process.warning); return; end
        end
        
        filesRawAll{iBlock,1} = filesRaw(:,1);
        filesRawAll{iBlock,2} = I;
    end
    
    %% Add OpenEphys to path
    addpath(parameters.path.ephys);
    
    %% STIMULUS ONSET DETECTION
    
    stimTimes = cell(nBlocks,2); % stimulus onset times
    
    % Different detection methods for stimulus onset times. This is heavily
    % dependent on the experimental paradigm.
    
    % User-supplied custom function
    if (~isempty_field(parameters,'parameters.general.stimPath'))
        cfg            = [];
        cfg.fpath      = parameters.general.stimPath; % Path to (custom) stimulus detection function
        cfg.loadpath   = loadPath;
        cfg.parameters = parameters;
        stimTimes = psr_wrapper_function(cfg);
    end
    
    %% Add FieldTrip to path
    if (parameters.process.lfp); [FT_FOUND,FT_PRESENT] = psr_ft_path(parameters,'add');
    else,                         FT_FOUND = false;
    end
    
    parameters.Fr = parameters.lfp.Fr;
    
    %% SPIKE + LFP DETECTION
    
    for iBlock = 1:nBlocks
        
        % Initialize
        filesRaw        = filesRawAll{iBlock,1};
        fileFlags       = filesRawAll{iBlock,2};
        loadPathBlock   = loadPath{iBlock};
        nFiles          = size(filesRaw,1);
        nProbes         = nFiles / parameters.general.nelectrodes;
        stimTimesTrial  = stimTimes(iBlock,:);
        stim            = parameters.general.stimuli{iBlock};
        ts_Spikes       = [];
        ts_LFP          = [];
        
        for iProbe = 1:nProbes
            
            % Check if TEMP file already exists
            if (~isempty(tempFilesFound) && ...
                    iProbe <= size(tempFilesFound,1) && ...
                    iBlock <= size(tempFilesFound,2))
                if (tempFilesFound(iProbe,iBlock)); continue; end
            end
            
            % Extract filenames for probe channels
            nChans     = parameters.general.nelectrodes;
            fileIdx    = (iProbe - 1) * nChans + (1:nChans);
            probeChans = filesRaw(fileIdx,1);
            keep       = fileFlags(fileIdx);
            chanIDs    = find(keep)';
            
            % Initialize
            artifactsProbe = [];
            dataChanRawSum = [];
            dataProbe      = [];
            timeProbe      = [];
            DATA_ERROR     = false;
            Fs_array       =  NaN(nChans,1);
            Ls_array       =  NaN(nChans,1);
            filenamesProbe = cell(nChans,1);
            
            for iChan = chanIDs
                
                % Load Open-Ephys data
                
                filename = [loadPathBlock probeChans{iChan,1}];
                filename = strtrim(filename);
                filenamesProbe{iChan} = filename;
                
                try % Load CONTINUOUS files [microvolts]
                    disp(['Loading ' filename '...']);
                    logInfo(logfilePath,'load_open_ephys_data_faster'); % Log the function call
                    [dataChanRaw, timestamps, info] = load_open_ephys_data_faster(filename);
                catch
                    logInfo(logfilePath,'load_open_ephys_data'); % Log the function call
                    [dataChanRaw, timestamps, info] = load_open_ephys_data(filename);
                end
                
                timestamps = timestamps - timestamps(1);
                
                if (parameters.develop.comparison) % Reduce data size for method comparison
                    if (~isempty_field(parameters,'parameters.develop.time'))
                        keep = 1:find(timestamps < parameters.develop.time,1,'last');
                        dataChanRaw = dataChanRaw(keep);
                        timestamps  = timestamps (keep);
                    end
                end
                
                dataChanRaw = single(dataChanRaw);
                
                % Check data properties
                
                Fs_array(iChan) = info.header.sampleRate;
                Ls_array(iChan) = length(dataChanRaw);
                [Fs,Ls] = checkDataProperties(Fs_array,Ls_array);
                if (isnan(Fs) || isnan(Ls))
                    str_1 = 'Inconsistent data found in file(s):';
                    str_2 = ['Directory: ' loadPathBlock];
                    psr_show_warning([str_1,probeChans',str_2],true);
                    DATA_ERROR = true;
                    break;
                end
                
                % Save sum of raw data
                
                if (parameters.ms.offset.run)
                    if (isempty(dataChanRawSum)); dataChanRawSum = dataChanRaw;
                    else,                         dataChanRawSum = dataChanRaw + dataChanRawSum;
                    end
                end
                
                % Magnetic stimulus artifact filter
                if (parameters.ms.denoise.raw.run)
                    logInfo(logfilePath,'psr_ms_denoise_raw'); % Log the function call
                    [dataChanRaw,artifacts] = psr_ms_denoise_raw(dataChanRaw,parameters,Fs);
                    artifactsProbe = [artifactsProbe;artifacts];
                end
                
                % Adaptive frequency spectrum filter
                if (parameters.filter.fft.run)
                    logInfo(logfilePath,'psr_artifact_fft'); % Log the function call
                    dataChanRaw = psr_artifact_fft(dataChanRaw,parameters,Fs);
                end
                
                %% SPIKING DATA
                
                if (parameters.process.spikes)
                    
                    % Band-pass filtering
                    cfg       = [];
                    cfg.Fs    = Fs;
                    cfg.order = parameters.spikes.bp.order;
                    cfg.lower = parameters.spikes.bp.lower;
                    cfg.upper = parameters.spikes.bp.upper;
                    
                    sLength   = length(dataChanRaw);
                    sPoints   = 60 * parameters.general.twin * Fs;
                    nSections = ceil(sLength / sPoints);
                    sPoints   = ceil(sLength / nSections);
                    itr       = 1;
                    
                    dataChanFiltered = zeros(1,sLength,'single');
                    for iSection = 1:nSections
                        I = itr:itr+sPoints-1; I(I > sLength) = [];
                        logInfo(logfilePath,'psr_bp_filter'); % Log the function call
                        dataChanFiltered(I) = psr_bp_filter(double(dataChanRaw(I)),cfg);
                        itr = I(end) + 1;
                    end
                    
                    if (isempty_field(ts_Spikes,'ts_Spikes.data')) % Initialize
                        ts_Spikes.data = zeros(nChans,length(dataChanFiltered),'int16');
                    end
                    
                    ts_Spikes.data(iChan,:) = int16(precision * dataChanFiltered);
                end
                
                %% LOCAL FIELD POTENTIAL
                
                if (parameters.process.lfp)
                    if (isempty(dataProbe) || isempty(timeProbe))
                        dataProbe = NaN(nChans,Ls);
                        timeProbe = NaN(nChans,Ls);
                    end
                    dataProbe(iChan,:) = dataChanRaw;
                    timeProbe(iChan,:) = timestamps;
                end
            end
            
            [Fs,Ls] = checkDataProperties(Fs_array,Ls_array);
            if (isnan(Fs) || isnan(Ls)); DATA_ERROR = true; end
            if (DATA_ERROR); continue; end
            
            % Data conversion and filtering
            
            if (parameters.process.lfp)
                
                parameters.Fs    = Fs;
                input.data       = dataProbe;
                input.timestamps = nanmean(timeProbe,1);
                logInfo(logfilePath,'psr_ft_convert2fieldtrip'); % Log the function call
                data = psr_ft_convert2fieldtrip(input,parameters);
                if (FT_FOUND)
                    logInfo(logfilePath,'psr_lfp_preprocessing'); % Log the function call
                    ts_LFP = psr_lfp_preprocessing(data,parameters);
                else
                    % TODO: Separate routine...
                end
            end
            
            if (parameters.process.spikes)
                ts_Spikes.Fs = Fs; % Hz, sampling rate of spike data
            end
            
            % Save
            
            metadata.filenames = filenamesProbe;
            metadata.artifacts = unique(artifactsProbe,'rows');
            metadata.subject   = subject;
            metadata.session   = session;
            metadata.stimtimes = stimTimesTrial;
            metadata.duration  = Ls / Fs;
            metadata.stimulus  = stim;
            metadata.probe     = iProbe;
            metadata.Fs        = Fs;
            
            if (parameters.ms.offset.run)
                logInfo(logfilePath,'psr_ms_detect_offset');
                metadata.stimoffset = psr_ms_detect_offset(dataChanRawSum,stimTimesTrial,Fs);
            end
            
            % Save temporary MAT file
            filename = [savePath ...
                'Temp_Probe_' num2str(iProbe,'%02d') ...
                '_Block_'     num2str(iBlock,'%02d') ...
                '.mat'];
            filesTemp{iProbe,iBlock} = filename;
            save(filename,'ts_Spikes','ts_LFP','metadata','parameters');
        end
    end
    
    % Remove FieldTrip from path
    if (FT_FOUND && ~FT_PRESENT); psr_ft_path(parameters,'remove'); end
    
    % Sort blocks
    
    nBlocks = size(filesTemp,2);
    
    sessions  = parameters.general.sessionIndex;
    nSessions = length(unique(sessions));
    trialIDs  = 1:nBlocks;
    itr       = 1;
    if (all(cellfun(@isnumeric,parameters.general.stimuli))) % Check if all stimuli are numbers
        for iSession = 1:nSessions
            stimSession = cell2mat(parameters.general.stimuli(sessions == iSession));
            [~,Isort] = sort(stimSession);
            n   = length(stimSession);
            I   = itr:itr+n-1;
            IDs = trialIDs(I);
            trialIDs(I) = IDs(Isort);
            itr = itr + n;
        end
    end
    filesTemp = filesTemp(:,trialIDs);
    
    % Log which files have been processed so far
    if (exist(logfilePath,'file')); save(logfilePath,'filesTemp','-append');
    else,                           save(logfilePath,'filesTemp');
    end
    
    % Remove OpenEphys from path
    rmpath(parameters.path.ephys);
    
end

%% Load temporary data processing file

filesSaved = [];
if (logfileFound || exist(logfilePath,'file'))
    [filesTemp,filesSaved] = psr_load_vars(logfilePath,{'filesTemp','filesSaved'});
end
if (isempty(filesSaved))
    filesSaved = cell(size(filesTemp,1),5); % Array for temporary data logging
end

% Fix array size
nProbes = size(filesTemp, 1);
mProbes = size(filesSaved,1);
d = nProbes - mProbes;
if     (d < 0); filesSaved = filesSaved(1:nProbes,:);
elseif (d > 0); filesSaved{nProbes,end} = [];
end

%% Force certain sections to be re-processed

nProbes    = size(filesSaved,1);
nSections  = size(filesSaved,2);
for iProbe = 1:nProbes
    for iSection = 1:nSections
        str = filesSaved(iProbe,iSection);
        I = any(strcmp(str,parameters.process.section));
        if (I || ~exist(filesSaved{iProbe,1},'file'))
            filesSaved{iProbe,iSection} = []; 
        end
    end
end

%% Merge metadata across trials

% Initialize arrays
nProbes = size(filesTemp,1);
nBlocks = size(filesTemp,2);

for iProbe = 1:nProbes
    
    if (strcmp(filesSaved{iProbe,2},'META')); continue; end
    
    metadata     = [];
    filenamesAll = cell(nBlocks,1);
    artifactsAll = cell(nBlocks,1);
    stimAllTypes = cell(nBlocks,1);
    stimAllAmps  =  NaN(nBlocks,1);
    onsetTimes   =  NaN(nBlocks,1);
    offsetTime   = 0;
    stimOffsets  = [];
    
    for iBlock = 1:nBlocks
        
        filename = filesTemp{iProbe,iBlock};
        
        if ~isempty(filename)
            
            metadata = psr_load_vars(filename,{'metadata'}); % Load temporary trial data
            
            if (~isempty(metadata))
                artifactsAll{iBlock} = metadata.artifacts;
                filenamesAll{iBlock} = metadata.filenames;
                
                if (iBlock == 1); stimAllTimes = cell(nBlocks,length(metadata.stimtimes)); end
                stimAllTimes(iBlock,:) = metadata.stimtimes;
                stimAllTypes{iBlock}   = metadata.stimulus;
                
                if (ischar(metadata.stimulus))
                    amp = regexp(metadata.stimulus,'\d*','Match');
                    amp = str2double(amp(:,1));
                else
                    amp = metadata.stimulus;
                end
                
                stimAllAmps(iBlock) = amp;
                
                if (~isempty_field(metadata,'metadata.stimoffset'))
                    stimOffsets(iBlock,:) = sort(metadata.stimoffset)';
                end
                
                onsetTime  = offsetTime;
                offsetTime = onsetTime + metadata.duration;
                onsetTimes(iBlock) = onsetTime;
            end
        end
    end
    
    if (~isempty(metadata))
        metadata.filenames  = filenamesAll;
        metadata.artifacts  = artifactsAll;
        metadata.stimtimes  = stimAllTimes;
        metadata.stimulus   = stimAllTypes;
        metadata.stimamps   = stimAllAmps;
        metadata.blockonset = onsetTimes;
        metadata.blockdur   = diff([onsetTimes;offsetTime]);
        metadata.duration   = offsetTime;
        if (~isempty(stimOffsets)); metadata.stimoffset = stimOffsets; end
        
        parameters.general.savelist = {'metadata'}; % What variables to save in output MAT file
        
        filesSaved{iProbe,1} = saveFile([],[],metadata,parameters,savePath,true); % Save to MAT file
    end
    
    filesSaved{iProbe,2} = 'META';
    save(logfilePath,'filesSaved','-append'); % Update output file array
end

% Combine stimulus offsets if needed
psr_ms_combine_offsets(filesSaved(:,1),parameters);

%% Local field potential

if (parameters.process.lfp)
    
    % Add FieldTrip to path
    if (parameters.process.lfp)
        [FT_FOUND,FT_PRESENT] = psr_ft_path(parameters,'add');
    end
    
    if (FT_FOUND)
        for iProbe = 1:nProbes
            
            if (strcmp(filesSaved{iProbe,3},'LFP')); continue; end
            disp(['LFP processing for probe ' num2str(iProbe,'%02d') '...']);
            
            % Initialize
            metadata  = [];
            dataArray = [];
            dataProbe = cell(0,0);
            
            % Load all trials for probe
            for iBlock = nBlocks:-1:1
                filename = filesTemp{iProbe,iBlock};
                if (~isempty(filename))
                    load(filename,'ts_LFP','metadata','parameters');
                    dataProbe{iBlock} = ts_LFP;
                end
            end
            
            if (~isempty(dataProbe))
                
                % Data filtering
                
                if (parameters.lfp.artifact.chan.run); dataProbe = psr_lfp_artifact_channel(dataProbe,parameters); end
                if (parameters.lfp.mean_subtract);     dataProbe = psr_lfp_mean_subtraction(dataProbe); end
                
                % Artifact removal
                
                artifacts = [];
                if (parameters.lfp.artifact.psd.run); artifacts.psd = psr_lfp_artifact_detection_psd(dataProbe,parameters); end
                if (parameters.lfp.artifact.amp.run); artifacts.amp = psr_lfp_artifact_detection_amp(dataProbe,parameters); end
                dataProbe = psr_lfp_artifact_removal(dataProbe,artifacts,parameters);
                
                % Convert LFP data to structure array
                
                for iBlock = nBlocks:-1:1
                    dataBlock = dataProbe{iBlock};
                    if (~isempty(dataBlock))
                        dataBlock.artifacts = dataProbe{iBlock}.artifacts;
                        fields = fieldnames(dataBlock);
                        for iField = 1:length(fields); dataArray(iBlock).(fields{iField}) = dataBlock.(fields{iField}); end
                    end
                end
                
            end
            
            % Save LFP output
            
            if (~isempty(metadata))
                parameters.general.savelist = {'freq'};
                saveFile([],dataArray,metadata,parameters,savePath,true);
            end
            filesSaved{iProbe,3} = 'LFP';
            save(logfilePath,'filesSaved','-append'); % Update output file array
        end
    end
    
    % Remove FieldTrip from path
    if (FT_FOUND && ~FT_PRESENT); psr_ft_path(parameters,'remove'); end
end

%% Spike processing

if (parameters.process.spikes)
    
    sortMethod = lower(parameters.sorting.method);
    
    %% Spike sorting
    
    for iProbe = 1:nProbes % Do clustering per probe across all stimulus conditions
        
        % Check if output file + temp file (after sorting) exist
        if (strcmp(filesSaved{iProbe,4},'SPK')); continue; end
        
        % Ignore missing blocks
        
        filesProbe = filesTemp(iProbe,:);
        keep = ~cellfun(@isempty,filesProbe);
        filesProbe = filesProbe(keep);
        nBlocks = length(filesProbe);
        
        % Check for missing temporary files
        
        MISSING_FILE = false;
        for iBlock = 1:nBlocks
            filename = filesProbe{iBlock};
            if (~exist(filename,'file')); MISSING_FILE = true; break; end
        end
        if (MISSING_FILE) % If not all files are present
            for iBlock = 1:nBlocks
                filename = filesProbe{iBlock};
                if (exist(filename,'file')); delete(filename); end % Delete temporary spikes file
            end
            disp(['Skipping probe ' num2str(iProbe) '...']);
            continue; % Move on to next probe
        end
        
        disp(['Spike sorting probe ' num2str(iProbe,'%02d') '...']);
        tSort = tic; % Measure duration of spike detection and sorting
        
        %% MERGE DATA ACROSS BLOCKS
        
        % Initialize
        dataProbe = cell(1,nBlocks);
        Fs_array  = zeros(nBlocks,1);
        Ls_array  = zeros(nBlocks,1);
        
        for iBlock = 1:nBlocks
            filename = filesProbe{iBlock};
            if ~isempty(filename)
                load(filename,'ts_Spikes');
                parameters.Fs     = ts_Spikes.Fs;
                Fs_array (iBlock) = ts_Spikes.Fs;
                Ls_array (iBlock) = size(ts_Spikes.data,2);
                dataProbe{iBlock} = ts_Spikes.data;
            end
        end
        
        dataProbe = cell2mat(dataProbe);
        
        %% Remove zero-signal channels
        
        del = all(dataProbe == 0,2);
        dataProbe(del,:) = [];
        nChans = size(dataProbe,1);
        
        %% Initialize and check data
        
        spikesAll = []; % structure to contain all spike information from all experimental blocks
        Fs = checkDataProperties(Fs_array); % Only check sampling frequency, since data length may vary across blocks
        
        if (nChans >= parameters.general.minchans && ~isnan(Fs))
            
            parameters.Fs = Fs;
            
            %% Calculate spike info
            
            spikesAll             = psr_sst_background_noise(spikesAll,dataProbe,parameters);
            spikesAll.Fs          = Fs;
            spikesAll.info.dur    = Ls_array ./ Fs_array;
            spikesAll.info.thresh = parameters.spikes.thresh * spikesAll.info.bgn;
            spikesAll.assigns     = [];
            spikesAll.spiketimes  = [];
            spikesAll.waveforms   = [];
            
            %% Spike detection
            
            if (any(strcmp(sortMethod,{'fmm','iso','kfm','ops','ost','spc','ums'}))) % Spike detection
                spikesAll = psr_sst_detection(spikesAll,dataProbe,parameters);
                spikesAll.info.detected = true;
            else
                spikesAll.info.detected = false;
            end
            
            %% SPIKE SORT
            
            switch sortMethod
                case 'cbp'; spikesAll = psr_sst_sorting_CBP(spikesAll,dataProbe,parameters);
                case 'fmm'; spikesAll = psr_sst_sorting_FMM(spikesAll,parameters);
                case 'iso'; spikesAll = psr_sst_sorting_ISO(spikesAll,parameters);
                case 'kfm'; spikesAll = psr_sst_sorting_KFM(spikesAll,parameters);
                case 'kst'; spikesAll = psr_sst_sorting_KST(spikesAll,dataProbe,parameters,savePath);
                case 'ops'; spikesAll = psr_sst_sorting_OPS(spikesAll,dataProbe,parameters);
                case 'ost'; spikesAll = psr_sst_sorting_OST(spikesAll,parameters);
                case 'spc'; spikesAll = psr_sst_sorting_SPC(spikesAll,parameters);
                case 'ums'; spikesAll = psr_sst_sorting_UMS(spikesAll,parameters);
            end
            
            if (parameters.develop.timing) % For development
                tEnd = toc(tSort); % Spike sorting time
                
                filePathTiming = [savePath
                    'timing_'    metadata.subject, ...
                    '_'          metadata.session, ...
                    '_P' num2str(metadata.probe,'%02d')];
                
                if (exist(filePathTiming,'file'))
                    load(filePathTiming);
                else
                    timings = [];
                    methods = [];
                end
                methods{end+1} = sortMethod;
                timings(end+1) = tEnd;
                save(filePathTiming,'methods','timings');
            end
            
            %% SET BLOCK INDEX FOR EACH SPIKE
            
            if (~isempty_field(spikesAll,'spikesAll.spiketimes'))
                
                % Initialize arrays
                blocks      = zeros(size(spikesAll.spiketimes),'int16');
                offsetTime  = 0;
                offsetSpike = 0;

                for iBlock = 1:nBlocks
                    load(filesProbe{iBlock},'metadata'); % Load temporary trial data
                    onsetTime   = offsetTime;
                    offsetTime  = onsetTime + metadata.duration;
                    onsetSpike  = offsetSpike + 1;
                    offsetSpike = find(spikesAll.spiketimes <= offsetTime,1,'last');
                    blocks(onsetSpike:offsetSpike) = iBlock;
                end

                spikesAll.blocks = blocks;
            end
        end
        
        %% SAVE
        
        metadata = [];
        if (~isempty(filesProbe))
            load(filesProbe{1},'metadata');  % Load metadata from first block
            parameters.general.savelist = {'spikes','parameters'}; % What variables to save in output MAT file
            saveFile(spikesAll,[],metadata,parameters,savePath,true); % Save to MAT file
        end
        filesSaved{iProbe,4} = 'SPK';
        save(logfilePath,'filesSaved','-append'); % Update output file array
        
    end
        
    %% Cluster processing
    
    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    for iProbe = 1:nProbes
        
        if (strcmp(filesSaved{iProbe,5},'CLS')); continue; end
        disp(['Processing clusters for probe ' num2str(iProbe,'%02d') '...']);
        
        % Initialize
        spikes = [];
        freq   = [];
        
        % Load spike data
        filename = filesSaved{iProbe,1};
        
        if (~isempty(filename))
            
            [spikes,spikes_old,metadata] = psr_load_vars(filename,{'spikes','spikes_old','metadata'});
            
            % Raw data files
            filesProbe = filesTemp(iProbe,:);
            keep = ~cellfun(@isempty,filesProbe);
            filesProbe = filesProbe(keep);
            
            if (~isempty_field(spikes,'spikes.spiketimes') && ~isempty(filesProbe))
                
                if (~isempty(spikes_old));       spikes = spikes_old; end
                if (~parameters.process.delete); spikes_old = spikes; end
                
                % Convert some variables
                spikes = psr_freq2spikes(spikes,freq);
                if (isfield(metadata,'artifacts')); spikes.info.artifacts.raw = metadata.artifacts;
                else,                               spikes.info.artifacts.raw = [];
                end
                
                spikes = psr_sst_assigns_reorder(spikes); % Re-order assign IDs
                spikes = psr_sst_cluster_quality(spikes,parameters);
                                
                if (parameters.ms.denoise.off.run);  spikes = psr_ms_denoise_off     (spikes,metadata,parameters);   end % Stimulus artifact removal
                if (parameters.filter.chan.rip.run); spikes = psr_sst_filter_chan_rip(spikes,parameters,filesProbe); end % Remove ripples on non-spike channels
                if (parameters.filter.chan.mse.run); spikes = psr_sst_filter_chan_mse(spikes,parameters);            end % Remove noise on non-spike channels
                if (parameters.filter.chan.loc.run); spikes = psr_sst_filter_chan_loc(spikes,parameters);            end
                
                spikes = psr_sst_cluster_quality   (spikes,parameters); % Pre-merge quality control
                spikes = psr_sst_cluster_remove    (spikes);            % Delete noise clusters
                spikes = psr_sst_spike_align       (spikes,parameters,filesProbe); % Centre on peak
                spikes = psr_sst_white_noise       (spikes,parameters); % Substitute tagged channels with white noise [again after aligning]
                spikes = psr_sst_features          (spikes,parameters); % Calculate low-dimensional features
                spikes = psr_sst_cluster_merge     (spikes,parameters); % Merge clusters
                spikes = psr_sst_cluster_quality   (spikes,parameters); % Post-merge quality control
                spikes = psr_sst_cluster_isolation (spikes,parameters); % Isolation quality of clusters
                spikes = psr_sst_cluster_thresholds(spikes,parameters); % Classify cluster quality
                spikes = psr_sst_filter_spikes     (spikes,parameters,'delete'); % Now filter using all enabled methods
                
                % Save
                spikes = orderfields(spikes);
                if (parameters.process.delete); save(filename,psr_varname(spikes),                        '-append');
                else,                           save(filename,psr_varname(spikes),psr_varname(spikes_old),'-append');
                end
            end
        end
        
        % Update output file array
        filesSaved{iProbe,5} = 'CLS';
        save(logfilePath,'filesSaved','-append');
    end
    warning('on', MSGID);
    
    disp(['Spike sorting completed. MAT file(s) saved to: "' savePath '"']);
    
end

%% Stability check (for development)
if (parameters.develop.comparison)
    filesSaved = psr_stability_check(filesSaved,filesTemp,sortMethod); %#ok
    save(logfilePath,'filesSaved','-append'); % Update output file array
end

%% Delete temporary files
if (parameters.process.delete)
    delete(logfilePath);
    for iProbe = 1:nProbes
        for iBlock = 1:nBlocks
            delete(filesTemp{iProbe,iBlock});
        end
    end
end

end

function filename = saveFile(spikes,freq,metadata,parameters,savePath,append)

if (nargin < 6); append = false; end

saveList = parameters.general.savelist;

filename = [savePath ...
    'PSR_' metadata.subject                  ...
    '_'    strjoin(metadata.session, '-')    ...
    '_P'   num2str(metadata.probe,   '%02d')];

if (~isempty(freq));   freq   = orderfields(freq);   end %#ok
if (~isempty(spikes)); spikes = orderfields(spikes); end %#ok
parameters = orderfields(parameters);                    %#ok
metadata   = orderfields(metadata);                      %#ok

filename = [filename '.mat'];
if (append && exist(filename,'file')); save(filename,saveList{:},'-append');
else,                                  save(filename,saveList{:});
end

end

function val = extractStringFromPath(filepath,token)
val      = 0;
filepath = fliplr(filepath);
k = strfind(filepath,'\');
k = k(2);
filepath = filepath(1:k);
filepath = fliplr(filepath);
elems = regexp(filepath, token, 'tokens', 'once');
if (~isempty(elems)); str = elems{1};
    if (~isempty(str)); val = str2double(elems{1}); end
end

end

function [Fs,Ls] = checkDataProperties(Fs_array,Ls_array)

if (nargin < 1); Fs_array = []; end
if (nargin < 2); Ls_array = []; end

Fs_array = Fs_array(~isnan(Fs_array));
Ls_array = Ls_array(~isnan(Ls_array));

Fs = NaN;
Ls = NaN;

if ~isempty(Fs_array) && all(Fs_array == Fs_array(1)); Fs = Fs_array(1); end
if ~isempty(Ls_array) && all(Ls_array == Ls_array(1)); Ls = Ls_array(1); end

end

function logInfo(logfilePath,str)

% Save info in logfile

callStack = cell(0,0);
if (exist(logfilePath,'file'))
    callStack = psr_load_vars(logfilePath,{'callStack'});
    callStack{end+1,1} = str;
    save(logfilePath,'callStack','-append');
else
    callStack{end+1,1} = str;
    save(logfilePath,'callStack');
end
end

%------------- END OF CODE --------------