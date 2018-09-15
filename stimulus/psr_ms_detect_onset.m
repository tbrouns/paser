function stimulusTimes = psr_ms_detect_onset(loadPath,parameters)

% PSR_MS_DETECT_ONSET - Detects timings of magnetic field stimulus
% This function finds the magnetic field stimulus onsets in the raw
% data, which are marked by peaks in the analog-to-digital
% converter (ADC) signal.
%
% Syntax:  stimulusTimes = psr_ms_detect_onset(loadPath,parameters)
%
% Inputs:
%    loadPath   - Path to ADC files [string]
%    parameters - See README and PSR_PARAMETERS_GENERAL
%
% Outputs:
%    stimulusTimes - Onset and offset of stimuli [sec]
%
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

nBlocks  = length(loadPath);
stimulusTimes = cell(nBlocks,1);

% Files to load
ext     = '.continuous';
pattern = 'ADC6';

for iBlock = 1:nBlocks
        
    % Load new file
    
    files = dir([loadPath{iBlock} '\*' pattern '*' ext]);
    files = char(files.name);
    
    if (size(files,1) == 1)
        file = files(1,:);
        file = strtrim(file);
    else
        continue;
    end
            
    % Filter raw data
                
    file = [loadPath{iBlock} file]; % Filename
    
    % Load CONTINUOUS files [microvolts]
    try    [signal, ~, info] = load_open_ephys_data_faster(file); 
    catch, [signal, ~, info] = load_open_ephys_data(file);
    end

    Fs = info.header.sampleRate; % Sampling rate in Hz
                
    %% Stimulus detection
    
    signal  = rescale(signal',-1,1);
    peaks   = find(signal > parameters.ms.detect.threshold); % check where magnetic pulse occurs
    offsets = find(diff(peaks) > 1); % offset of pulse
    onsets  = offsets + 1; % take onset of pulse
    onsets  = [peaks(1),       peaks(onsets)]; % first index always an onset
    offsets = [peaks(offsets), peaks(end)]; % last index always an offset
    onsets  = ( onsets - 1) / Fs;
    offsets = (offsets - 1) / Fs;
    
    % Threshold pulse duration
    
    dur = offsets - onsets;
    id  = dur > parameters.ms.detect.min_dur;
    onsets  =  onsets(id)';
    offsets = offsets(id)';
    
    %% Display results

    period = mean(diff(onsets));
    disp(['Detected ' num2str(length(onsets)) ' magnetic stimulus onsets at ' num2str(1 / period) ' Hz.']);
    
    %% Save
    
    stimulusTimes{iBlock} = [onsets,offsets];
    
end

end

%------------- END OF CODE --------------