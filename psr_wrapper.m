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
%    parameters - See PSR_PARAMETERS_GENERAL and PSR_BATCH_PROCESSING
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

%% Load processing parameters
parameters = psr_load_parameters(parameters);

%% Constants
nChans    =    parameters.general.nelectrodes;
pattern   =    parameters.general.rawpattern;
ext       =    parameters.general.extension;
precision = 10^parameters.general.precision;
tempFilePath = [savePath 'Temp_vars.mat'];

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

fileTemp = dir(tempFilePath);
fileTemp = char(fileTemp.name);
tempFileFound = ~isempty(fileTemp);

%% Start data processing

if (~tempFileFound || any(~tempFilesFound(:)))
    
    %% Extract experimental conditions from folder name
    
    parameters.general.stimuli = cell(nBlocks,1);
    
    for iBlock = 1:nBlocks
        filepath = loadPath{iBlock};
        if (~psr_isempty_field(parameters,'parameters.general.blockpattern'))
            parameters.general.stimuli{iBlock} = extractStringFromPath(filepath,['_(\d*)' parameters.general.blockpattern]);
        else
            string = regexp(filepath,'\','split');
            parameters.general.stimuli{iBlock} = string{end-1};
        end
    end
    
    %% Find raw data files
    
    filesUnsorted = cell(nBlocks,1);
    for iBlock = 1:nBlocks
        files = dir([loadPath{iBlock} '\*' pattern '*' ext]);
        if (size(files,1) == 0); return; end
        filesUnsorted{iBlock} = char(files.name);
    end
    
    %% Check raw data files and sort them in correct chronological order
    
    filesRawAll = cell(nBlocks,1);
    
    for iBlock = 1:nBlocks
        nFiles = length(filesUnsorted{iBlock}(:,1));
        filesRaw = cell(nFiles,2);
        for iFile = 1:nFiles
            filename   = filesUnsorted{iBlock}(iFile,:);
            filename   = strtrim(filename);
            itr        = strfind(filename,pattern) + length(pattern);
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
        
        filename = filesRaw{1,2}; % Compare first file with others
        I = strcmp(filename,filesRaw(:,2));
        
        if (any(~I))
            strMain = ['ERROR: raw data filenames are not consistent. Check data folder: "' loadPath{iBlock} '"'];
            str = repmat('-',1,length(strMain));
            disp(str); 
            disp(strMain);
            disp(str)
            return; 
        end
        
        filesRawAll{iBlock} = filesRaw; % All data files have similar data format
    end
    
    %% Add OpenEphys to path
    addpath(parameters.path.ephys);
    
    %% STIMULUS ONSET DETECTION
    
    stimTimes = cell(nBlocks,2); % stimulus onset times
    
    % Different detection methods for stimulus onset times. This is heavily
    % dependent on the experimental paradigm.
       
    % User-supplied custom function
    if (~psr_isempty_field(parameters,'parameters.general.stimPath'))
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
        filesRaw        = filesRawAll{iBlock};
        nFiles          = length(filesRaw);
        nProbes         = nFiles / parameters.general.nelectrodes;
        stimTimesTrial  = stimTimes(iBlock,:);
        stim            = parameters.general.stimuli{iBlock};
        ts_Spikes       = [];
        ts_LFP          = [];
        
        for iProbe = 1:nProbes
            
            % Initialize
            artifactsProbe = [];
            
            % Check if TEMP file already exists
            if (~isempty(tempFilesFound) && ...
                    iProbe <= size(tempFilesFound,1) && ...
                    iBlock <= size(tempFilesFound,2))
                if (tempFilesFound(iProbe,iBlock)); continue; end
            end
            
            for iChan = 1:nChans
                
                % Load Open-Ephys data
                
                iFile = (iProbe - 1) * nChans + iChan;
                
                filename = [loadPath{iBlock} filesRaw{iFile}];
                filename = strtrim(filename);
                
                try % Load CONTINUOUS files [microvolts]
                    disp(['Loading ' filename '...']);
                    [dataChanRaw, timestamps, info] = load_open_ephys_data_faster(filename);
                    timestamps = timestamps / info.header.sampleRate;
                catch
                    [dataChanRaw, timestamps, info] = load_open_ephys_data(filename);
                end
                
                dataChanRaw = single(dataChanRaw);
                
                % Store variables
                if (iChan == 1) % Initialize
                    Fs_array = zeros(nChans,1);
                    Ls_array = zeros(nChans,1);
                end
                
                Fs = info.header.sampleRate;
                Fs_array(iChan) = Fs; % Sampling frequencies should be equal
                Ls_array(iChan) = length(dataChanRaw);
                
                if (parameters.ms.offset.run)
                    if (iChan == 1); dataChanRawSum = dataChanRaw;
                    else,            dataChanRawSum = dataChanRaw + dataChanRawSum;
                    end
                end
                
                % Magnetic stimulus artifact filter
                if (parameters.ms.denoise.raw.run)
                    [dataChanRaw,artifacts] = psr_ms_denoise_raw(dataChanRaw,parameters,Fs);
                    artifactsProbe = [artifactsProbe;artifacts];
                end
                
                % Adaptive frequency spectrum filter
                if (parameters.filter.fft.process); dataChanRaw = psr_artifact_fft(dataChanRaw,parameters,Fs); end
                
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
                    
                    dataChanFiltered = zeros(1,sLength);
                    for iSection = 1:nSections
                        I = itr:itr+sPoints-1; I(I > sLength) = [];
                        dataChanFiltered(I) = psr_bp_filter(double(dataChanRaw(I)),cfg);
                        itr = I(end) + 1;
                    end
                    
                    if (iChan == 1); ts_Spikes.data = zeros(nChans,length(dataChanFiltered),'int16'); end
                    
                    ts_Spikes.data(iChan,:) = int16(precision * dataChanFiltered);
                end
                
                %% LOCAL FIELD POTENTIAL
                
                if (parameters.process.lfp)
                    if (iChan == 1)
                        dataProbe = zeros(nChans,Ls_array(iChan));
                        timeProbe = zeros(nChans,Ls_array(iChan));
                    end
                    dataProbe(iChan,:) = dataChanRaw;
                    timeProbe(iChan,:) = timestamps - timestamps(1);
                end
                
            end
            
            % Check data
            
            [Fs,Ls] = checkDataProperties(Fs_array,Ls_array);
            if (isempty(Fs) || isempty(Ls)); continue; end
            
            % Data conversion and filtering
            
            if (parameters.process.lfp)
                parameters.Fs    = Fs;
                input.data       = dataProbe;
                input.timestamps = mean(timeProbe);
                data = psr_ft_convert2fieldtrip(input,parameters);
                if (FT_FOUND)
                    ts_LFP = psr_lfp_preprocessing(data,parameters);
                else
                    % TODO: Separate routine...
                end
            end
            
            if (parameters.process.spikes)
                ts_Spikes.Fs = Fs; % Hz, sampling rate of spike data
            end
            
            % Save
                        
            metadata.artifacts = sortrows(artifactsProbe);
            metadata.subject   = subject;
            metadata.session   = session;
            metadata.stimtimes = stimTimesTrial;
            metadata.duration  = Ls / Fs;
            metadata.stimulus  = stim;
            metadata.probe     = iProbe;
            metadata.Fs        = Fs;
            
            if (parameters.ms.offset.run)
                metadata.stimoffset = psr_ms_detect_offset(dataChanRawSum,stimTimesTrial,Fs);
            end
            
            % Save temporary MAT file
            filename = [savePath ...
                'Temp_Probe_' num2str(iProbe,  '%02d') ...
                '_Block_'     num2str(iBlock,  '%02d') ...
                '.mat'];
            filesTemp{iProbe,iBlock} = filename;
            save(filename,'ts_Spikes','ts_LFP','metadata','parameters');
        end
    end
    
    % Remove FieldTrip from path
    if (FT_FOUND && ~FT_PRESENT); psr_ft_path(parameters,'remove'); end
    
    % Temporarily save all workspace variables. Useful when continuing after early stopping
    
    keepvars = {...
        'parameters',    ...
        'filesTemp',     ...
        'metadata',      ...
        'loadPath',      ...
        'savePath',      ...
        'tempFilePath'};
    clearvars('-except', keepvars{:});
    save(tempFilePath);
else
    
    % Load temporary data processing file
    filesSaved = [];
    load(tempFilePath);
    
    % Force certain sections to be re-processed
    parameters = psr_load_parameters(parameters); % TEMP
    nProbes    = size(filesSaved,1);
    nSections  = size(filesSaved,2);
    if (~isempty(filesSaved) && ~psr_isempty_field(parameters,'parameters.process.section'))
        for iProbe = 1:nProbes
            for iSection = 1:nSections
                str = filesSaved(iProbe,iSection);
                I = any(strcmp(str,parameters.process.section));
                if (I); filesSaved{iProbe,iSection} = []; end
            end
        end
    end
    
end

%% Remove OpenEphys from path
rmpath(parameters.path.ephys);

%% Sort trials

nProbes = size(filesTemp,1);
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

if (~exist('filesSaved','var')) % Initialize array if not loaded
    filesSaved = cell(nProbes,5); % Output files
end

%% Merge metadata across trials

% Initialize arrays
nProbes = size(filesTemp,1);
nBlocks = size(filesTemp,2);

for iProbe = 1:nProbes
    
    if (exist(filesSaved{iProbe,1},'file') && strcmp(filesSaved{iProbe,2},'META')); continue; end
    
    artifactsAll =  cell(nBlocks,1);
    stimAllTypes =  cell(nBlocks,1);
    stimAllAmps  = zeros(nBlocks,1);
    onsetTimes   = zeros(nBlocks,1);
    offsetTime   = 0;
    stimOffsets  = [];
    
    for iBlock = 1:nBlocks
        
        load(filesTemp{iProbe,iBlock},'metadata','parameters'); % Load temporary trial data
        
        artifactsAll{iBlock} = metadata.artifacts;
        
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
        
        if (isfield(metadata,'stimoffset')); stimOffsets(iBlock,:) = sort(metadata.stimoffset)'; end
        
        onsetTime  = offsetTime;
        offsetTime = onsetTime + metadata.duration;
        onsetTimes(iBlock) = onsetTime;
        
    end
    
    metadata.artifacts  = artifactsAll;
    metadata.stimtimes  = stimAllTimes;
    metadata.stimulus   = stimAllTypes;
    metadata.stimamps   = stimAllAmps;
    metadata.trialonset = onsetTimes;
    metadata.duration   = offsetTime;
    if (~isempty(stimOffsets)); metadata.stimoffset = stimOffsets; end
    
    parameters.general.savelist = {'metadata'}; % What variables to save in output MAT file
    
    filesSaved{iProbe,1} = saveFile([],[],metadata,parameters,savePath,true); % Save to MAT file
    filesSaved{iProbe,2} = 'META';
    save(tempFilePath,'filesSaved','-append'); % Update output file array
            
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
            
            if (exist(filesSaved{iProbe,1},'file') && strcmp(filesSaved{iProbe,3},'LFP')); continue; end
            
            disp(['LFP processing for probe ' num2str(iProbe,'%02d') '...']);
            
            % Load all trials for probe
            dataProbe = cell(nBlocks,1);
            for iBlock = 1:nBlocks
                load(filesTemp{iProbe,iBlock},'ts_LFP','metadata');
                dataProbe{iBlock} = ts_LFP;
            end
            
            load(filesTemp{iProbe,1},'parameters'); % Load parameters from first trial
            
            parameters = psr_load_parameters(parameters); % TEMP
            
            % Data filtering
            
            if (parameters.lfp.artifact.chan.run); dataProbe = psr_lfp_artifact_channel(dataProbe,parameters); end
            if (parameters.lfp.mean_subtract);     dataProbe = psr_lfp_mean_subtraction(dataProbe); end
                
            % Artifact removal
            
            artifacts = [];
            if (parameters.lfp.artifact.psd.run); artifacts.psd = psr_lfp_artifact_detection_psd(dataProbe,parameters); end
            if (parameters.lfp.artifact.amp.run); artifacts.amp = psr_lfp_artifact_detection_amp(dataProbe,parameters); end
            dataProbe = psr_lfp_artifact_removal(dataProbe,artifacts,parameters);
            % TEMP: NEEDS TESTING
            
            % Convert LFP data to structure array
            
            freqArray = [];
            for iBlock = nBlocks:-1:1
                freq           = dataProbe{iBlock};
                freq.artifacts = dataProbe{iBlock}.artifacts;
                fields = fieldnames(freq);
                for iField = 1:length(fields); freqArray(iBlock).(fields{iField}) = freq.(fields{iField}); end
            end
            
            % Save LFP output
            
            parameters.general.savelist = {'freq'};
            filesSaved{iProbe,1} = saveFile([],freqArray,metadata,parameters,savePath,true);
            filesSaved{iProbe,3} = 'LFP';
            save(tempFilePath,'filesSaved','-append'); % Update output file array
        end
    end
    
    % Remove FieldTrip from path
    if (FT_FOUND && ~FT_PRESENT); psr_ft_path(parameters,'remove'); end
end

%% Spike sorting

parameters = psr_load_parameters(parameters);

if (parameters.process.spikes)
    
    sortMethod = lower(parameters.sorting.method);
    
    for iProbe = 1:nProbes % Do clustering per probe across all stimulus conditions
        
        % Check if output file + temp file (after sorting) exist
        if (exist(filesSaved{iProbe,1},'file') && strcmp(filesSaved{iProbe,4},'SPK')); continue; end
        
        % Check for missing temporary files
        
        MISSING_FILE = false;
        for iBlock = 1:nBlocks
            if (~exist(filesTemp{iProbe,iBlock},'file')); MISSING_FILE = true; break; end
        end
        if (MISSING_FILE) % If not all files are present
            for iBlock = 1:nBlocks
                filename = filesTemp{iProbe,iBlock};
                if (exist(filename,'file')); delete(filename); end % Delete temporary spikes file
            end
            disp(['Skipping probe ' num2str(iProbe) '...']);
            continue; % Move on to next probe
        end
        
        disp(['Spike sorting probe ' num2str(iProbe,'%02d') '...']);
        tSort = tic; % Measure duration of spike detection and sorting
        
        %% MERGE DATA ACROSS TRIALS
        
        % Initialize
        dataProbe = [];
        Fs_array  = zeros(nBlocks,1);
        Ls_array  = zeros(nBlocks,1);
        
        for iBlock = 1:nBlocks
            
            load(filesTemp{iProbe,iBlock});
            
            parameters.Fs    = ts_Spikes.Fs;
            Fs_array(iBlock) = ts_Spikes.Fs;
            Ls_array(iBlock) = size(ts_Spikes.data,2);
            
            if (parameters.develop.timing) % For development
                parameters = psr_load_parameters(parameters);
                n = round(parameters.Fs * parameters.develop.time);
                ts_Spikes.data = ts_Spikes.data(:,1:n);
            end
            
            dataProbe = [dataProbe,ts_Spikes.data]; %#ok
        end
                
        %% Calculate spike info
        
        Fs = checkDataProperties(Fs_array);
        if (isempty(Fs)); continue; end
        parameters.Fs = Fs;
        
        spikesAll             = []; % structure to contain all spike information from all trials
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
            filename = [savePath 'timings.mat'];
            if (exist(filename,'file'))
                load(filename);
            else
                timings = [];
                methods = [];
            end
            methods{end+1} = sortMethod;
            timings(end+1) = tEnd;
            save([savePath 'timings.mat'],'methods','timings');
            continue;
        end
        
        %% SET BLOCK INDEX FOR EACH SPIKE
                
        % Initialize arrays
        
        nBlocks     = size(filesTemp,2);
        blocks      = zeros(size(spikesAll.spiketimes),'int16');
        offsetTime  = 0;
        offsetSpike = 0;
        
        for iBlock = 1:nBlocks
            load(filesTemp{iProbe,iBlock},'metadata','parameters'); % Load temporary trial data
            onsetTime   = offsetTime;
            offsetTime  = onsetTime + metadata.duration;
            onsetSpike  = offsetSpike + 1;
            offsetSpike = find(spikesAll.spiketimes <= offsetTime,1,'last');
            blocks(onsetSpike:offsetSpike) = iBlock;
        end
        
        spikesAll.blocks = blocks;
        
        %% SAVE
                
        load(filesTemp{iProbe,1},'metadata'); % Load metadata from first trial
        parameters.general.savelist = {'spikes','parameters'}; % What variables to save in output MAT file
        filesSaved{iProbe,1} = saveFile(spikesAll,[],metadata,parameters,savePath,true); % Save to MAT file
        filesSaved{iProbe,4} = 'SPK';
        save(tempFilePath,'filesSaved','-append'); % Update output file array
        
    end
    
    if (parameters.develop.timing); return; end    
        
    %% Cluster processing
    
    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    for iProbe = 1:nProbes
        
        if (exist(filesSaved{iProbe,1},'file') && strcmp(filesSaved{iProbe,5},'CLS')); continue; end
        disp(['Processing clusters for probe ' num2str(iProbe,'%02d') '...']);
        
        spikes = [];
        freq   = [];
        
        filename = filesSaved{iProbe,1};
        load(filename);
        
        parameters.probeID = iProbe; % TEMP (for ripple visualization)
        
        if (psr_isempty_field(spikes,'spikes.spiketimes')); continue; end
        
        filesProbe = filesTemp(iProbe,:);
        
        % Convert some variables
        spikes = psr_sst_freq2spikes(spikes,freq);
        if (isfield(metadata,'artifacts')); spikes.info.artifacts.raw = metadata.artifacts;
        else,                               spikes.info.artifacts.raw = []; 
        end
        
        % Noise removal
        
        parameters = psr_load_parameters(parameters);
        
        if (parameters.ms.denoise.off.run);       spikes = psr_ms_denoise_off    (spikes,metadata,parameters);   end
        if (parameters.filter.spikes.ripple.run); spikes = psr_sst_remove_ripples(spikes,parameters,filesProbe); end % Remove ripples on non-spike channels
        if (parameters.filter.spikes.noise.run);  spikes = psr_sst_remove_noise  (spikes,parameters);            end % Remove noise on non-spike channels
                
        spikes = psr_sst_cluster_quality   (spikes,parameters); % Pre-merge quality control
        spikes = psr_sst_cluster_remove    (spikes);            % Delete noise clusters
        spikes = psr_sst_spike_align       (spikes,parameters,filesProbe); % Centre on peak
        spikes = psr_sst_features          (spikes,parameters); % Calculate low-dimensional features
        spikes = psr_sst_cluster_merge     (spikes,parameters); % Merge clusters
        spikes = psr_sst_cluster_quality   (spikes,parameters); % Post-merge quality control 
        spikes = psr_sst_cluster_isolation (spikes,parameters); % Isolation quality of clusters    
        spikes = psr_sst_cluster_thresholds(spikes,parameters); % Classify cluster quality
        spikes = psr_sst_filter_spikes     (spikes,parameters,'delete'); % Now delete RPVs, if enabled
        
        % Save 
        spikes = orderfields(spikes);
        save(filename,psr_varname(spikes),'-append');
        
        % Update output file array
        filesSaved{iProbe,5} = 'CLS';
        save(tempFilePath,'filesSaved','-append'); 
    end
    warning('on', MSGID);
    
    disp(['Spike sorting completed. MAT file(s) saved to: "' savePath '"']);
    
end

%% Stability check (for development)
if (parameters.develop.comparison)
    filesSaved = psr_stability_check(filesSaved,filesTemp,sortMethod); %#ok
    save(tempFilePath,'filesSaved','-append'); % Update output file array
    return;
end

%% Delete temporary files
if (parameters.process.delete)
    delete(tempFilePath);
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
if (append && exist(filename,'file') > 0); save(filename,saveList{:},'-append');
else,                                      save(filename,saveList{:});
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

Fs = [];
Ls = [];

if ~isempty(Fs_array)
    if (any(Fs_array ~= Fs_array(1)))
        disp('Sampling frequency mismatch');
    else
        Fs = Fs_array(1);
    end
end

if ~isempty(Ls_array)
    if (any(Ls_array ~= Ls_array(1)))
        disp('Data length mismatch');
    else
        Ls = Ls_array(1);
    end
end
end

%------------- END OF CODE --------------