function [FT_FOUND,FT_PRESENT] = psr_ft_path(parameters,method)

% PSR_FT_PATH - Adds or removes the FieldTrip toolbox to search path
%
% Syntax:  [FT_FOUND,FT_PRESENT] = psr_ft_path(parameters,method)
%
% Inputs:
%    parameters - See README
%    method     - String which can be set to:
%                 'add'    :    Adds the FieldTrip toolbox to path
%                 'remove' : Removes the FieldTrip toolbox to path
%
% Outputs:
%    FT_FOUND   - Boolean indicating whether the FieldTrip toolbox has been
%                 found and added to search path
%    FT_PRESENT - Boolean indicating whether the FieldTrip is already on
%                 the search path and didn't have to be added
%
% See also: PSR_WRAPPER, FT_DEFAULTS

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% Initialize
FT_FOUND = true;

if (strcmp(method,'add')) % Add FieldTrip path
        
    % Check if FieldTrip (FT) exists on path
    names = who('global');
    k     = [strcmp(names,'ft_default')]; %#ok
    if ~isempty(find(k,1)); FT_PRESENT = true; % FieldTrip already exists on path
    else,                   FT_PRESENT = false;
    end
    
    if (~FT_PRESENT)
        if (isempty_field(parameters,'parameters.path.ft'))
            FT_FOUND = false; % FT not found and no path specified
            disp('FieldTrip path has not been set.');
        else
            % Add FieldTrip to path
            MSGID1 = 'MATLAB:dispatcher:pathWarning'; warning('off',MSGID1); 
            MSGID2 = 'MATLAB:rmpath:DirNotFound';     warning('off',MSGID2);
            addpath(parameters.path.ft);
            try
                ft_defaults
                [~, ftpath] = ft_version;
                rmpath(fullfile(ftpath, 'external', 'signal'));
                rmpath(fullfile(ftpath, 'external', 'stats'));
                rmpath(fullfile(ftpath, 'external', 'image'));
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