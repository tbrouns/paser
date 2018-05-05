function metadata = psr_cam_sync(metadata,parameters)

% Find offset between NIDAQ and ADC

Fs         = parameters.Fs;
Fr         = 0.5 * 1000 / parameters.analysis.sync.bin;
halfWin    = round(0.5 * Fs * parameters.analysis.sync.bin    / 1000);
desyncMax  = round(      Fs * parameters.analysis.sync.desync / 1000);
chanID     = parameters.analysis.sync.chan;
nSessions  = length(metadata.session);
sessionIDs = metadata.sessionIndex;

metadata.camoffsets = cell(nSessions,1);

for iSession = 1:nSessions
    
    loadPath = [parameters.general.loadPath '\' metadata.session{iSession}];
    loadPathSub = [loadPath '\NIDAQ\'];
    
    filesNIDAQ = dir([loadPathSub 'Data_*.mat']);
    filesNIDAQ = char(filesNIDAQ.name);
    nfiles     = size(filesNIDAQ,1);
    
    if (nfiles == 0); continue; end
    
    % Combine NIDAQ data
    
    NIDAQ.analog = cell(nfiles,1);
    NIDAQ.time   = cell(nfiles,1);
    
    for iFile = 1:nfiles
        
        filename = filesNIDAQ(iFile,:);
        filepath = [loadPathSub filename];
        load(filepath,'Data');
        
        I = regexp(filename, 'Data_(\d*).', 'tokens', 'once');
        I = str2double(cell2mat(I));
        
        NIDAQ.analog{I} = Data.Analog(:,chanID);
        NIDAQ.time  {I} = Data.Time;
        
    end

    NIDAQ.analog = cat(1, NIDAQ.analog{:});
    NIDAQ.time   = cat(1, NIDAQ.time  {:});

    [NIDAQ_onset,NIDAQ_offset] = psr_transient_detection(NIDAQ.analog,parameters.analysis.sync.thresh);
    
    NIDAQ_onset  = round(Fs * NIDAQ.time(NIDAQ_onset))  + 1; 
    NIDAQ_offset = round(Fs * NIDAQ.time(NIDAQ_offset)) + 1; 
    
    % Remove faulty timings
    
    del_on  = diff(NIDAQ_onset)  < 0;
    del_off = diff(NIDAQ_offset) < 0;
    del     = find(del_on | del_off) + 1;
    NIDAQ_onset (del) = [];
    NIDAQ_offset(del) = [];
            
    % Get ADC times
    
    blockIDs = (metadata.sessionIndex == iSession);
    stimTypes = metadata.stimtimes(blockIDs,2);
    I = strcmp(stimTypes,'interval');
    stimTimes = metadata.stimtimes(I,1);
    stimTimes = cat(1, stimTimes{:});
    
    ADC_onset  = round(Fs * stimTimes(:,1)) + 1;
    ADC_offset = round(Fs * stimTimes(:,2)) + 1;
    
    % Find approximate offset between signals
    
    nLength = round(Fs * metadata.duration);
    
    signalNIDAQ = zeros(nLength,1);
    signalADC   = zeros(nLength,1);
    
    window = -halfWin:halfWin;
    
    NIDAQ_win_on  = bsxfun(@plus,NIDAQ_onset, window);
    NIDAQ_win_off = bsxfun(@plus,NIDAQ_offset,window);
    
    ADC_win_on  = bsxfun(@plus,ADC_onset, window);
    ADC_win_off = bsxfun(@plus,ADC_offset,window);
    
    signalNIDAQ(NIDAQ_win_on)  =  1;
    signalNIDAQ(NIDAQ_win_off) = -1;
    
    signalADC(ADC_win_on)  =  1;
    signalADC(ADC_win_off) = -1;
        
    % Remove initial zeros
    
    iniOffsetNIDAQ = find(signalNIDAQ,1,'first');
    iniOffsetADC   = find(signalADC,  1,'first');
    
    signalDiff = iniOffsetNIDAQ - iniOffsetADC;
            
    % Get exact offset for start of each video recording
    
    d = ADC_onset + signalDiff - NIDAQ_onset';
    [~,I] = min(abs(d),[],2);
    
    I = sub2ind(size(d),1:length(ADC_onset),I');
    d = d(I);
    d(abs(d) > desyncMax) = NaN;
    metadata.camoffsets{iSession} = (signalDiff - d)' / Fs;
                   
end

end