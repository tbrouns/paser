function [FT_FOUND,FT_PRESENT] = psr_ft_path(parameters,method)

% Initialize
FT_FOUND = true;

if (strcmp(method,'add')) % Add field trip path
        
    % Check if FieldTrip (FT) exists on path
    names = who('global');
    k     = [strcmp(names,'ft_default')]; %#ok
    if ~isempty(find(k,1)); FT_PRESENT = true; % FieldTrip already exists on path
    else,                   FT_PRESENT = false;
    end
    
    if (~FT_PRESENT)
        if (psr_isempty_field(parameters,'parameters.path.ft'))
            FT_FOUND = false; % FT not found and no path specified
            disp('FieldTrip path has not been set.');
        else
            % Add FieldTrip to path
            MSGID1 = 'MATLAB:dispatcher:pathWarning'; warning('off',MSGID1); 
            MSGID2 = 'MATLAB:rmpath:DirNotFound';     warning('off',MSGID2);
            addpath(parameters.path.ft);
            try
                ft_defaults
                rmpath(fullfile(parameters.path.ft, 'external', 'signal'))
                rmpath(fullfile(parameters.path.ft, 'external', 'stats'))
                rmpath(fullfile(parameters.path.ft, 'external', 'image'))
            catch
                FT_FOUND = false; % Path set, but FT not found
                disp('FieldTrip not found on given path.');
            end
            warning('on',MSGID1);
            warning('on',MSGID2);
        end
    end
elseif (strcmp(method,'remove'))
    try % Remove FT from path
        psr_remove_path(parameters.path.ft);
        clearvars -global
    catch
        disp('Unable to remove FieldTrip from path.');
    end
end

end