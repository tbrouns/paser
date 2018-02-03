function tf = psr_exist_in_file(filename,varname) 
tf = false;
if (exist(filename,'file'))
    variableInfo = who('-file', filename);
    tf = ismember(varname, variableInfo);
end
end