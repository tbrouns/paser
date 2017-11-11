function MFAtimes = psr_mfa_detection(loadPath,parameters)

nTrials  = length(loadPath);
MFAtimes = cell(nTrials,2);

% Files to load
ext         = '.continuous';
pattern     = 'ADC';
fileIndices = [4 7];
FILE_FLAG   = true;
nFiles = length(fileIndices);

for iTrial = 1:nTrials
    
    % Load new files
    files = cell(nFiles,1);
    
    for iFile = 1:nFiles
        
        file = dir([loadPath{iTrial} '\*' pattern num2str(fileIndices(iFile)) '*' ext]);
        file = char(file.name);
        
        if (size(file,1) == 1)
            file = file(1,:);
            file = strtrim(file);
        else
            FILE_FLAG = false;
            break;
        end
        
        files{iFile} = file;
    end
        
    if (~FILE_FLAG); continue; end % next trial
    
    % Initialize arrays
    
    data_channels =  cell(nFiles,1);
    Fs_array      = zeros(nFiles,1);
    
    % Filter raw data
    
    for iFile = 1:nFiles
        
        file = [loadPath{iTrial} files{iFile}]; % Filename
        
        try % Load CONTINUOUS files [microvolts]
            [data_ADC, ~, info] = load_open_ephys_data_faster(file);
        catch
            [data_ADC, ~, info] = load_open_ephys_data(file);
        end
        
        Fs = info.header.sampleRate; % Sampling rate in Hz
        Fs_array(iFile) = Fs; % Save in array, because we have multiple sampling frequencies (should all be equal)
        
        % Band-pass filter
        
        [B,A] = butter(parameters.mfa.bp_order,[parameters.mfa.bp_lower parameters.mfa.bp_upper]/(Fs/2),'bandpass');
        data_ADC = filtfilt(B,A,data_ADC); % Zero-phase digital filtering
        
        % Normalize
        
        data_ADC = data_ADC / max(data_ADC);
        data_channels{iFile} = single(data_ADC); % convert to single type
    end
    
    parameters.Fs = mode(Fs_array); % Just take most frequently occuring sampling frequency
    
    if (parameters.general.stims(iTrial) > 0) % If we have a stimulus
        [MFAtimes{iTrial,1},MFAtimes{iTrial,2}] = mfa_detection(data_channels,parameters);
    else % Randomly select MFA times in order to get non-stimulus control data
        dur = (length(data_channels{1}) - 1) / Fs;
        MFAtimes{iTrial,1} = sort(dur * rand(1,parameters.mfa.control))';
        MFAtimes{iTrial,2} = sort(dur * rand(1,parameters.mfa.control))';
    end
end

end


function [spikeTimes_1,spikeTimes_2] = mfa_detection(signals,parameters)

% PSR_MFA_DETECTION - Detects location of magnetic field artifacts (MFAs)
% in signal, which are marked by large peaks in the analog-to-digital
% converter (ADC) signal. Typically, there are two intertwined MFA signals,
% so we end up with two different output MFA arrays. The first of these is
% used as the stimulus onset time.
%
% Syntax:  [spikeTimes_1,spikeTimes_2] = psr_mfa_detection(signals,parameters)
%
% Inputs:
%    signals    - Cell array, in which each cell contains a one-dimensional ADC time series
%    parameters - See "PSR_MFA_DETECTION" section in PSR_PARAMETER_DEFAULT
%
% Outputs:
%    spikeTimes_1 - Locations of MFA peaks of 1st intertwined signal [sec]
%    spikeTimes_2 - Locations of MFA peaks of 2nd intertwined signal [sec]
%
% See also: PSR_WRAPPER

% PASER: Processing and Analysis Schemes for Extracellular Recordings
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept.
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

Fs = parameters.Fs;
window_size = (parameters.mfa.p2ptime / 1000) * Fs; % Maximum number of samples between min-max peaks

nChans = length(signals);
spikeTimesAll = cell(nChans,1);

for iChan = 1:nChans
    
    signalChannel = signals{iChan};
    signalChannel = signalChannel';
    
    % Calculate threshold
    stdev     = psr_mad(signalChannel); % median absolute deviation
    threshold = parameters.mfa.thresh * stdev;
    
    % Find peaks and troughs
    [pks_min,loc_min] = findpeaks(double(-signalChannel)); % Locations of local minima
    [pks_max,loc_max] = findpeaks(double( signalChannel)); % Locations of local maxima
    pks_min           = -pks_min;
    
    % Calculate p2p amplitudes of adjacent min-max pair
    
    Nmin = size(pks_min,2);
    Nmax = size(pks_max,2);
    
    if (Nmin > Nmax)
        pks_min = pks_min(1:Nmax);
        loc_min = loc_min(1:Nmax);
    elseif (Nmax > Nmin)
        pks_max = pks_max(1:Nmin);
        loc_max = loc_max(1:Nmin);
    end
    
    if (loc_min(1) < loc_max(1))
        loc_1 = loc_min;
        loc_2 = loc_max;
        pks_1 = pks_min;
        pks_2 = pks_max;
    else
        loc_1 = loc_max;
        loc_2 = loc_min;
        pks_1 = pks_max;
        pks_2 = pks_min;
    end
    
    % Detect putative artifacts based on amplitudes
    
    % Each local minimum or maximum has one or two neighbours, so we need
    % to iterate over the various pair-wise combinations
    
    X = [2  1  1];
    x = [0  0 -1];
    Y = fliplr(X);
    y = fliplr(x);
    
    loc_all = [];
    
    for i = 1:length(X)
        
        dn  =   abs(loc_1(X(i):end+x(i))-loc_2(Y(i):end+y(i)));    % Distance between min and max
        loc = mean([loc_1(X(i):end+x(i));loc_2(Y(i):end+y(i))],1); % Between-peak locations
        
        p1  = pks_1(X(i):end+x(i));
        p2  = pks_2(Y(i):end+y(i));
        
        % Determine which are the minima and which are the maxima
        if (sum(p1) < sum(p2)); pmin = p1; pmax = p2;
        else,                   pmin = p2; pmax = p1;
        end
        
        % Filter based on p2p time
        id   = (dn <= window_size);
        pmax = pmax(id);
        pmin = pmin(id);
        loc  =  loc(id);
        
        % Filter based on max and min amplitudes
        id      = pmax > threshold & pmin < -threshold; % Signal is symmetric
        loc_all = [loc_all, loc(id)]; %#ok, artifact locations
        
    end
    
    % Peak locations in seconds
    spikeTimes = (loc_all - 1) / Fs;
    spikeTimes = sort(spikeTimes);
    
    periods = diff(spikeTimes);
    periodSorted = sort(periods);
    I = floor(0.5 * length(periods));
    
    if (I > 1)
        periodLower = median(periodSorted(1:I));
        periodUpper = median(periodSorted(I:end));
        [~,periodIDs] = min([abs(periods - periodLower);abs(periods - periodUpper)]);
        spikeTimesAll{iChan} = spikeTimes(periodIDs == 1); % Extract dual spikes and only take first
    else
        spikeTimesAll{iChan} = [];
    end
end

%% Check where ADC signals have peak overlap

[X,Y] = meshgrid(spikeTimesAll{2},spikeTimesAll{1});
[D,I] = min(abs(X-Y));
I = I(D < (window_size / Fs));
id    = false(length(spikeTimesAll{1}),1);
id(I) = true; % Part of main signal

spikeTimes_1 = spikeTimesAll{1}( id)'; % If overlap, then part of 1st signal
spikeTimes_2 = spikeTimesAll{1}(~id)'; % Otherwise, part of 2nd signal

%% Display results

period_1 = mean(diff(spikeTimes_1));
period_2 = mean(diff(spikeTimes_2));

N1 = floor(((length(signalChannel) - 1) / Fs) * (1 / period_1));
N2 = floor(((length(signalChannel) - 1) / Fs) * (1 / period_2));

disp(['Signal 1: detected ' num2str(length(spikeTimes_1)) ' of ' num2str(N1) ' magnetic field artifacts at ' num2str(1 / period_1) ' Hz.']);
disp(['Signal 2: detected ' num2str(length(spikeTimes_2)) ' of ' num2str(N2) ' magnetic field artifacts at ' num2str(1 / period_2) ' Hz.']);

%------------- END OF CODE --------------

end