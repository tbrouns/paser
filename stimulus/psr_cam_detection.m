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
    
    if (size(files,1) > 1)
        continue;
    else
        file = files(1,:);
        file = strtrim(file);
    end
        
    % Filter raw data
    pathToFile = [loadPath{iTrial} file]; % Filename

    try % Load CONTINUOUS files [microvolts]
        [signal, ~, info] = load_open_ephys_data_faster(pathToFile);
    catch
        [signal, ~, info] = load_open_ephys_data(pathToFile);
    end

    Fs = info.header.sampleRate; % Sampling rate in Hz

    % Normalize

    signal = signal / max(signal);
    signal = single(signal); % convert to single type
                
    %% Camera onset detection
        
    signalDiff = diff(signal);    
    
    % On- and offsets
    transitions = find(abs(signalDiff) > threshold);
   
    % Join neighbouring transitions
    d = diff(transitions);
    transitions = transitions(find(d ~= 1) + 1);
    
    % Check if we are dealing with onset or offset
    amplitudes = signalDiff(transitions);
    amplitudes = amplitudes ./ abs(amplitudes);
        
    id = find(diff(amplitudes) == -2);
    onset  = transitions(id);
    offset = transitions(id + 1);
    
    % Save times
    camTimes{iTrial} = [onset,offset] / Fs;
    
end

end