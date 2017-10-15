function signalFiltered = psr_artifact_fft(signal,parameters,Fs)

% PSR_ARTIFACT_FFT - Filters raw time series signal in its power spectrum.
% This function detects high intensity peaks in the signal's power spectrum
% using median absolute deviation, subsequently removes them and returns
% the filtered signal.
%
% Syntax:  [dataFiltered] = psr_artifact_fft(data,parameters,Fs)
%
% Inputs:
%    signal     - Vector of time series signal
%    parameters - See "PSR_ARTIFACT_FFT" section in PSR_PARAMETER_DEFAULT
%    Fs         - Sampling frequency of data [Hz]
%
% Outputs:
%    signalFiltered - Filtered time series signal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: PSR_WRAPPER, PSR_PARAMETER_DEFAULT

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% email address: t.s.n.brouns@gmail.com
% Date: 2017

%------------- BEGIN CODE --------------

nlength = length(signal); % Number of data points
fwin    = round((parameters.filter.fft_freq / Fs) * nlength); % Moving frequency window

Y      = fft(signal); % Perform fast-fourier transform
Y_abs  = abs(Y); % Only take real absolute part of power spectrum to do thresholding
fstart = 1; % Frequency starting point for FFT window
peaks  = []; % Location of peaks in power spectrum
flag   = true; % condition to determine 

while flag
   
    fend = fstart + fwin; % Update frequency window
    
    if (fend > nlength) % Check if end of signal is reached
        fstart = nlength - fwin; 
        fend   = nlength;
        flag   = false; % Last window to process
    end
    
    Y_win = Y_abs(fstart:fend); % Extract FFT window
    thresh = parameters.filter.fft_thresh * psr_mad(Y_win); % Calculate threshold
    
    del = Y_win < thresh; % Do thresholding
    I   = find(~del); % Keep raw locations of thresholded data
    Y_win(del) = []; % Ignore sub-threshold data
    
    if (length(Y_win) > 3) % Criterion is needed for 'findpeaks'
        [~,loc] = findpeaks(Y_win); % Detect all peaks
        loc = I(loc); % Find locations of peaks in raw data
    else
        loc = I; % Just take all locations
    end
        
    peaks = [peaks;loc+fstart-1]; %#ok  
    
    fstart = fstart + round(0.5 * fwin); % Update starting frequency
        
end

fwin = ceil((parameters.filter.fft_pad / Fs) * nlength); % Size of window to remove artifact

peaks = unique(peaks); 
peaks = sort(peaks);
peaks = bsxfun(@plus,peaks,-fwin:fwin); % Create windows to remove artifacts
peaks = peaks(:);
peaks(peaks < 1)             = [];
peaks(peaks > length(Y_abs)) = [];
Y(peaks) = 0; % Remove artifacts from original power spectrum

signalFiltered = ifft(Y,'symmetric'); % Inverse FFT to acquire filtered signal

%------------- END OF CODE --------------

end