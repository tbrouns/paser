function camTimes = psr_cam_detection(loadPath)

nTrials   = length(loadPath);
camTimes  = cell(nTrials,1);
threshold = 0.1; % One-tenth of maximum derivative

% Files to load
ext     = '.continuous';
pattern = 'ADC3';

for iTrial = 1:nTrials
    
    %% Load new files
    
    files = dir([loadPath{iTrial} '\*' pattern '*' ext]);
    files = char(files.name);
    
    if (size(files,1) == 1)
        file = files(1,:);
        file = strtrim(file);
    else
        continue;
    end
        
    % Filter raw data
    pathToFile = [loadPath{iTrial} file]; % Filename

    try % Load CONTINUOUS files [microvolts]
        [signal, ~, info] = load_open_ephys_data_faster(pathToFile);
    catch
        [signal, ~, info] = load_open_ephys_data(pathToFile);
    end

    Fs = info.header.sampleRate; % Sampling rate in Hz

    [onset,offset] = psr_transient_detection(signal,threshold);
    
    % Save times
    camTimes{iTrial} = [onset,offset] / Fs;
    
end

end