function [data,nTrials] = psr_lfp_conversion(data)

% "Data" can be given as single time-series (Nchans x Npoints), or as a
% cell array of such matrices or FieldTrip data structures

if (iscell(data))
    nTrials = length(data);
else % Single time-series (Nchans x Npoints)
    data = {data};
    nTrials = 1;
end

% Check for FieldTrip, and do conversion if necessary

for iTrial = 1:nTrials
    if isfield(data{iTrial},'trial') 
        data{iTrial} = data{iTrial}.trial{1};
    end
end

end