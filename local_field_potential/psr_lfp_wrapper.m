function output = psr_lfp_wrapper(inputs,parameters)

% Initialize LFP structure
output = [];

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
else
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
    
    if (strcmp(inputs.method,'tfa')) % time-frequency analysis
        
        freq       = []; %#ok
        data       = inputs.data;
        timestamps = inputs.timestamps;
        stimtimes  = inputs.stimtimes;
        method     = stimtimes{2};
        
        data = psr_convert2fieldtrip(data,timestamps,stimtimes,parameters);
        [data,parameters] = psr_lfp_preprocessing(data,parameters); % FT preprocessing
        if (~isempty(data) && strcmp(method,'onset')); freq = psr_lfp_timefreq(data,parameters);
        else,                                          freq = data;
        end
        
        output.freq = freq;
        output.parameters = parameters;
    
    elseif (strcmp(inputs.method,'plot'))
    
        freq = inputs.freq;
        cfg  = inputs.cfg;
        psr_lfp_plotting(freq,cfg);
        
    end
    
    if (~FT_PRESENT)
        try % Remove FT from path
            psr_remove_path(parameters.path.ft);
            clearvars -global
        catch
            disp('Unable to remove FieldTrip from path.');
        end
    end
    
    disp('LFP processing completed.');
end

end
