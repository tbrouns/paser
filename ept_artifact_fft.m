function dataFiltered = ept_artifact_fft(data,parameters,Fs)

nlength = length(data);
fwin = round((parameters.filter.fft_freq / Fs) * nlength);

Y      = fft(data);
Y_abs  = abs(Y);
fstart = 1;
freq   = [];
flag   = true;

while flag
   
    fend = fstart + fwin;
    
    if (fend > nlength)
        fstart = nlength - fwin; 
        fend   = nlength;
        flag = false; 
    end
    
    Y_win = Y_abs(fstart:fend);
    thresh = parameters.filter.fft_thresh * ept_mad(Y_win);
    
    del = Y_win < thresh;
    I   = find(~del); % keep raw locations of thresholded data
    Y_win(del) = []; % ignore sub-threshold data
    
    if (length(Y_win) > 3) % criterion is needed for 'findpeaks'
        [~,loc] = findpeaks(Y_win); 
        loc = I(loc); % find locations of peaks in raw data
    else
        loc = I;
    end
        
    freq = [freq;loc+fstart-1]; %#ok  
    
    fstart = fstart + round(0.5 * fwin);
        
end

fwin = ceil((parameters.filter.fft_pad / Fs) * nlength);

freq = unique(freq);
freq = sort(freq);
freq = bsxfun(@plus,freq,-fwin:fwin);
freq = freq(:);
freq(freq < 1)         = [];
freq(freq > length(Y_abs)) = [];
Y(freq) = 0;

dataFiltered = ifft(Y,'symmetric');

%% PLOTTING
% figure;plot(1:length(Y_abs),Y_abs); hold on; scatter(freq,zeros(size(freq)));
% Y_abs(freq) = 0;
% figure;plot(1:length(Y_abs),Y_abs); hold on; scatter(freq,zeros(size(freq)));
% 
% figure; plot(1:length(data),data);
% figure; plot(1:length(dataFiltered),dataFiltered);

end