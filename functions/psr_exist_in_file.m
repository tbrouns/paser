function tf = psr_exist_in_file(filepath,varname) 

% PSR_EXIST_IN_FILE - Check if variable exists in MAT file
%
% Syntax:  tf = psr_exist_in_file(filepath,varname) 
%
% Inputs:
%    filename - Path to MAT file
%    varname  - Name of variable of interest
%
% Outputs:
%    tf - Boolean which will equal true if variable exists in MAT file and
%         will equal false otherwise
% 
% Example:
%   var = 'This is a variable'; % Create a dummy variable
%   save('file.mat','var'); % Save the dummy
%   tf = psr_exist_in_file('file.mat','var'); % Check if the dummy exists
%   disp(tf); % Yes, it exists in the file!
% 
% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

% 

tf = false;
if (exist(filepath,'file'))
    variableInfo = who('-file', filepath);
    tf = ismember(varname, variableInfo);
end
end