function varargout = psr_load_vars(filepath,varnames)

    % filepath: location of file you want to load
    % varname: name of variable that needs to be loaded from file

    MSGID = 'MATLAB:load:variableNotFound';
    warning('off', MSGID);
    output = load(filepath,varnames{:});
    nvars = length(varnames);
    varargout = cell(nvars,1);
    for iVar = 1:nvars
        varname = varnames{iVar};
        if (psr_isfield(output,['output.' varname]))
            varargout{iVar} = output.(varname);
        else
            varargout{iVar} = [];
        end
    end
    warning('on', MSGID);
        
end