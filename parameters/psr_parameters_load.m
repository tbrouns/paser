function parameters = psr_parameters_load(parameters,type)

% PSR_PARAMETERS_LOAD - Load parameters from default or custom script 
%
% Syntax:  parameters = psr_parameters_load(parameters,type)
%
% Inputs:
%    parameters - Should contain a field: "parameters.(type).configPath",
%                 where "type" is one of the three input types given below.
% 
%                 "configPath" should point to your custom parameter script
%                 for that particular type. See README. 
% 
%    type       - Can be 'general', 'analysis' or 'display'
%
% Outputs:
%    parameters - Structure containing the recently loaded parameters
%
% See also: PSR_PARAMETERS_GENERAL, PSR_PARAMETERS_ANALYSIS, PSR_PARAMETERS_DISPLAY

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (nargin < 1); parameters = [];  end
if (nargin < 2); type = 'general'; end
fname = ['parameters.' type '.configPath'];

tf = false;
if (isempty_field(parameters,fname))
    fprintf(['Path to config file for ' type ' parameters not set. '])
else
    try    
        run(parameters.(type).configPath);
        disp(['Loaded ' type ' parameters from "' parameters.(type).configPath '"']);
        tf = true;
    catch ME
        str = ['Error found in config file for ' type ' parameters:'];
        psr_show_warning({str,ME.message},true);
    end
end

if (~tf)
    fprintf(['Check "' fname '" field. Using default parameters...\n']);
    switch type
        case 'general';  psr_parameters_general; 
        case 'analysis'; psr_parameters_analysis;
        case 'display';  psr_parameters_display;
    end
end

parameters = orderfields(parameters);

end