function data = psr_ft_nan_removal(data)

% Deal with NaNs in data

nTrials = size(data.trial,2);
for iTrial = 1:nTrials
    Y = data.trial{iTrial};
    T = data.time {iTrial};
    nChans = size(Y,1);
    for iChan = 1:nChans
        
        % First sweep: Interpolate NaNs
        t = T;
        y = Y(iChan,:);
        del = isnan(y);
        y(del) = [];
        t(del) = [];
        Y(iChan,:) = interp1(t,y,T);    
        
        % Second sweep: Set any remaining NaNs to zero
        y = Y(iChan,:);
        del = isnan(y);
        Y(iChan,del) = 0;
        
    end
    data.trial{iTrial} = Y;
end

end