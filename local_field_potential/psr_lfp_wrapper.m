function [freq,parameters] = psr_lfp_wrapper(data,timestamps,stimtimes,parameters)

if exist('ft_default','var'); FT_PRESENT = true;
else,                         FT_PRESENT = false;
end

if (~FT_PRESENT)
    % Add FieldTrip to path
    addpath(parameters.path.ft);
    ft_defaults
end

data = psr_convert2fieldtrip(data,timestamps,stimtimes,parameters);
[data,parameters] = psr_preprocessing(data,parameters); % FT preprocessing
if (~isempty(data)); freq = psr_timefreq_analysis(data,parameters); end

if (~FT_PRESENT)
    % Remove FieldTrip path
    fpath = path;
    fpath = strsplit(fpath,';');
    k = ~cellfun(@isempty,strfind(fpath,parameters.path.ft));
    fpath = fpath(1,k);
    fpath = strjoin(fpath,';');
    rmpath(fpath);
    clearvars -global
end

end
