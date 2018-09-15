function data = psr_ft_nan_removal(data)

% PSR_FT_NAN_REMOVAL - Remove NaNs in data by interpolation
% Eliminate NaNs in the FieldTrip structure for LFP data. We first try to
% interpolate between start and end points, and otherwise substitute with
% zeros.
%
% Syntax:  data = psr_ft_nan_removal(data)
%
% Inputs:
%    data - A FieldTrip structure for LFP data
%
% Outputs:
%    data - Data without NaNs, and adds an "missing" field that indicates
%           which data points used to be NaNs
%
% See also: PSR_FT_CONVERT2FIELDTRIP, PSR_LFP_PREPROCESSING, PSR_LFP_TFA

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (isempty(data)); return; end

% Deal with NaNs in data

nTrials = size(data.trial,2);
for iTrial = 1:nTrials
    Y = data.trial{iTrial};
    T = data.time {iTrial};
    nChans = size(Y,1);
    M = false(nChans,length(T)); % missing data
    
    for iChan = 1:nChans
        
        % First sweep: Interpolate NaNs
        t = T;
        y = Y(iChan,:);
        del = isnan(y);
        y(del) = [];
        t(del) = []; 
        M(iChan,del) = true; % Tag missing data points
        
        % Interpolation requires at least two sample points in each dimension
        if (length(y) > 1 && length(t) > 1)
            Y(iChan,:) = interp1(t,y,T);   
        end
        
        % Second sweep: Set any remaining NaNs to zero
        y = Y(iChan,:);
        del = isnan(y);
        M(iChan,del) = true; % Tag missing data points
        Y(iChan,del) = 0;
        
    end
    data.trial  {iTrial} = Y;
    data.missing{iTrial} = M;
end

end