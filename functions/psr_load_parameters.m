function parameters = psr_load_parameters(parameters,type)

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