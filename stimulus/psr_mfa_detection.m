function MFAtimes = psr_mfa_detection(loadPath,parameters)

% PSR_MFA_DETECTION - Detects onset of magnetic field artifacts.
% This function finds the magnetic field artifact (MFA) onsets in the raw
% data, which are marked by large peaks in the analog-to-digital
% converter (ADC) signal.
%
% Syntax:  MFAtimes = psr_mfa_detection(loadPath,parameters)
%
% Inputs:
%    loadPath   - Path to ADC files [string]
%    parameters - See "PSR_MFA_DETECTION" section in PSR_PARAMETER_DEFAULT
%
% Outputs:
%    MFAtimes - Onset and offset of MFAs [sec]
%
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

nTrials  = length(loadPath);
MFAtimes = cell(nTrials,1);

% Files to load
ext         = '.continuous';
pattern     = 'ADC6';

for iTrial = 1:nTrials
        
    % Load new file
    
    files = dir([loadPath{iTrial} '\*' pattern '*' ext]);
    files = char(files.name);
    
    if (size(files,1) == 1)
        file = files(1,:);
        file = strtrim(file);
    else
        continue;
    end
            
    % Filter raw data
                
    file = [loadPath{iTrial} file]; % Filename

    try % Load CONTINUOUS files [microvolts]
        [signal, ~, info] = load_open_ephys_data_faster(file);
    catch
        [signal, ~, info] = load_open_ephys_data(file);
    end

    Fs = info.header.sampleRate; % Sampling rate in Hz
                
    %% MFA detection
    
    signal  = rescale(signal',-1,1);
    peaks   = find(signal > parameters.mfa.threshold); % check where magnetic pulse occurs
    offsets = find(diff(peaks) > 1); % offset of pulse
    onsets  = offsets + 1; % take onset of pulse
    onsets  = [peaks(1),       peaks(onsets)]; % first index always an onset
    offsets = [peaks(offsets), peaks(end)]; % last index always an offset
    onsets  = ( onsets - 1) / Fs;
    offsets = (offsets - 1) / Fs;
    
    % Threshold pulse duration
    
    dur = offsets - onsets;
    id  = dur > parameters.mfa.min_dur;
    onsets  =  onsets(id)';
    offsets = offsets(id)';
    
    %% Display results

    period = mean(diff(onsets));
    N = floor(((length(signal) - 1) / Fs) * (1 / period));
    disp(['Detected ' num2str(length(onsets)) ' of ' num2str(N) ' magnetic field artifacts at ' num2str(1 / period) ' Hz.']);
    
    %% Save
    
    MFAtimes{iTrial} = [onsets,offsets];
    
%     else % Randomly select MFA times in order to get non-stimulus control data
%         dur = (length(data_channels{1}) - 1) / Fs;
%         MFAtimes{iTrial,1} = sort(dur * rand(1,parameters.mfa.control))';
%         MFAtimes{iTrial,2} = sort(dur * rand(1,parameters.mfa.control))';
%     end
    
end

end

%------------- END OF CODE --------------