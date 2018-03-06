function tf = psr_exist_in_file(filename,varname) 

% Check if variable exists in MAT file

tf = false;
if (exist(filename,'file'))
    variableInfo = who('-file', filename);
    tf = ismember(varname, variableInfo);
end
end