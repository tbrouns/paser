function output = psr_load_var(filepath,varname)

    % filepath: location of file you want to load
    % varname: name of variable that needs to be loaded from file

    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    output = load(filepath,varname);
    if (psr_isfield(output,['output.' varname]))
        output = output.(varname);
    else
        output = [];
    end
    warning('on', MSGID);
        
end