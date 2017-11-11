function [freq,parameters] = psr_lfp_wrapper(data,timestamps,stimtimes,parameters)

% Initialize LFP structure
freq = [];

% Check FieldTrip path 

FT_PRESENT = false; % Check if FieldTrip (FT) exists on path
SKIP_LFP   = false;

if (isempty(parameters.path.ft))
    names = who('global');
    k     = [strcmp(names,'ft_default')]; %#ok
    if ~isempty(find(k, 1)); FT_PRESENT = true; % FT found, so don't add
    else
        SKIP_LFP = true; % FT not found and no path set
        disp('FieldTrip path not set.');
    end
elseif (~FT_PRESENT)
    % Add FieldTrip to path
    MSGID = 'MATLAB:dispatcher:pathWarning';
    warning('off',MSGID);
    addpath(parameters.path.ft);
    try 
        ft_defaults
    catch
        SKIP_LFP = true; % Path set, but FT not found
        disp('FieldTrip not found on given path.');
    end
    warning('on',MSGID);
end

% Do LFP processing

if (SKIP_LFP)
    disp('Skipping LFP processing...');
else
    data = psr_convert2fieldtrip(data,timestamps,stimtimes,parameters);
    [data,parameters] = psr_lfp_preprocessing(data,parameters); % FT preprocessing
    if (~isempty(data)); freq = psr_lfp_timefreq(data,parameters); end
    
    if (~FT_PRESENT)
        try % Remove FT from path
            fpath = path;
            fpath = strsplit(fpath,';');
            k = ~cellfun(@isempty,strfind(fpath,parameters.path.ft));
            fpath = fpath(1,k);
            fpath = strjoin(fpath,';');
            rmpath(fpath);
            clearvars -global
        catch
            disp('Unable to remove FieldTrip from path.');
        end
    end
    
    disp('LFP processing completed.');
end

end
