function [filenames,stims] = ept_load_files_session(filenames,loadPath)

% Save filenames in cell array where position depends on tetrode number and
% stimulus amplitude

if (nargin < 1); disp('No input files'); return; end
if (nargin < 2); loadPath = []; end
    
if (~isempty(loadPath) && loadPath(end) ~= '\'); loadPath = [loadPath, '\']; end

stims = cell(0,0);
numfiles = size(filenames,1);
filenamesCell = cell(0,0);
for iFile = 1:numfiles
    filename = strtrim(filenames(iFile,:));
    load([loadPath filename],'metadata');
    i = metadata.tetrode;
    j = [];
    for iStims = 1:length(stims)
        if (min(ismember(stims{iStims}, metadata.stimulus)))
           j = iStims; 
           break;
        end
    end
    if (isempty(j))
        filenamesCell{i,end+1} = filename; %#ok
        stims{end+1} = metadata.stimulus;  %#ok
    else
        filenamesCell{i,j} = filename;
    end
    
end
filenames = filenamesCell;
stims = cell2mat(stims);

end