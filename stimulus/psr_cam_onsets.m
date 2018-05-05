function metadata = psr_cam_onsets(metadata,parameters)

% loadPath: where to load annotated data from

% Check data
if (~psr_isfield(metadata,'metadata.camoffsets')); return; end

nSessions  = length(parameters.general.session);
touchTimes = cell(nSessions,1);

metadata.touchtimes = cell(nSessions,1);

for iSession = 1:nSessions
    
    if (isempty(metadata.camoffsets{iSession})); continue; end
    
    % Extract event times from ADC
    
    stimTypes = metadata.stimtimes(:,2);
    I = strcmp(stimTypes,'interval');
    stimTimes = metadata.stimtimes(I,1);
    stimTimes = cat(1, stimTimes{:});
    stimTimes = stimTimes + metadata.camoffsets{iSession};
         
    % Extract event data from NIDAQ
    
    foldername  =  parameters.general.session{iSession};
    loadPathRaw = [parameters.general.loadPath '\' foldername];
    
    loadPathEvents = [loadPathRaw '\General.mat'];
    
    if (~exist(loadPathEvents,'file')); continue; end
    
    load(loadPathEvents);
    eventTimes = {CGSave.Events.NI.Events.Time};
    eventTimes = cat(1, eventTimes{:});
    
    % Extract PointGrey data
    
    loadPathSub  = [loadPathRaw '\PointGrey\'];
    filesCamTemp = dir([loadPathSub 'Data_*.mat']);
    filesCamTemp = char(filesCamTemp.name);
    nFiles       = size(filesCamTemp,1);
    filesCam     = cell(nFiles,1);
    
    for iFile = 1:nFiles
        filename = filesCamTemp(iFile,:);
        filepath = [loadPathSub filename];
        index = regexp(filename, 'Data_(\d*).mat', 'tokens', 'once');
        index = cell2mat(index);
        if ~isempty(index)
            index = str2double(index);
            filesCam{index} = strtrim(filepath);
        end
    end
           
    % Extract clicked data
    
    loadPathSub = [parameters.general.loadPathClicked '\' foldername '\'];
    filesClickedTemp = dir([loadPathSub '\Data_*_Annotations.mat' ]);
    filesClickedTemp = char(filesClickedTemp.name);
    nFiles           = size(filesClickedTemp,1);

    touchTimesSession = cell(0,0);
    
    for iFile = 1:nFiles
        
        % Extract touch times
        
        filename = filesClickedTemp(iFile,:);
        filepath = [loadPathSub filename];
        index = regexp(filename, 'Data_(\d*)_', 'tokens', 'once');
        index = str2double(cell2mat(index));
        
        touchTimesSession{index,2} = strtrim(filename); % Save filename
        
        CurvesByFrame = [];
        load(filepath);
        idx = find(~cellfun(@isempty,CurvesByFrame));     
        nFrames = length(idx);
        types = false(nFrames,1);
        whiskerIDs = cell(0,0);
        for iFrame = 1:nFrames 
            frameInfo = CurvesByFrame{idx(iFrame)};
            nClicks = length(frameInfo);
            for iClick = 1:nClicks
                clickInfo = frameInfo{iClick};
                whiskerID = [clickInfo{5} num2str(clickInfo{6})];
                type = clickInfo{4};
                I = find(strcmp(whiskerIDs,whiskerID));
                if (isempty(I))
                    I = length(whiskerIDs) + 1;
                    whiskerIDs{I} = whiskerID;
                end
                if strcmp(type,'touch'); types(iFrame,I) = true; end
            end
        end
        
        nWhiskers = size(types,2);
        touchTimesRel = cell(nWhiskers,1);
        for iWhisker = 1:nWhiskers
            I = types(:,iWhisker);
            touchTimesRel{iWhisker} = idx(I);
        end
                
        % Link with PointGrey data
        
        Data = [];
        nFrames = length(CurvesByFrame);
        
        % Load corresponding PointGrey file
    
        if length(filesCam) >= index; fileCam = filesCam{index};
        else,                         fileCam = [];
        end
        
        if isempty(fileCam); continue; 
        else,                load(fileCam);
        end        
        
        if (psr_isempty_field(Data,'Data.Time')); continue; end
        if (nFrames ~= size(Data.Time,1))
            warning('Frame number mismatch');
        end
        
        videoStartTimeAbs = Data.Time(1,2); % Master clock
        [~,I] = min(abs(eventTimes(:,3) - videoStartTimeAbs)); % Link (NIDAQ) event times with (PointGrey) video start times

        videoStartTimeRel = eventTimes(I,1); % Relative to start of session

        [d,I] = min(abs(stimTimes(:,1) - videoStartTimeRel));
        thresh = parameters.analysis.sync.desync / 1000;
        
        if (d <= thresh)
            videoRange = stimTimes(I,:); % Start and end time of video recording
            for iWhisker = 1:nWhiskers
                touchRangeWhisker = touchTimesRel{iWhisker}; % Start and end time of touch
                % Find touch ranges
                if (~isempty(touchRangeWhisker))
                    touchTimeDiff = diff(touchRangeWhisker);
                    time_off = find(touchTimeDiff > 1);
                    time_on  = time_off + 1;
                    time_off = [time_off;length(touchRangeWhisker)];
                    time_on  = [1;time_on];
                    time_off = touchRangeWhisker(time_off);
                    time_on  = touchRangeWhisker(time_on);

                    % Touch times relative to start of recording
                    time_off = Data.Time(time_off,1) - Data.Time(1,1); 
                    time_on  = Data.Time(time_on, 1) - Data.Time(1,1); 

                    touchRangeWhisker = [time_on,time_off];
                    
                    dt = diff(videoRange);
                    if (all(time_on <= dt))
                        touchTimesSession{index}{iWhisker,1} = videoRange;
                        touchTimesSession{index}{iWhisker,2} = touchRangeWhisker;
                        touchTimesSession{index}{iWhisker,3} = whiskerIDs{iWhisker};
                    end
                end
            end
        end
    end
    
    I = cellfun(@isempty,touchTimesSession(:,1));
    touchTimesSession(I,:) = [];
    
    touchTimes{iSession} = touchTimesSession;
    
end

metadata.touchtimes = touchTimes;
    
end