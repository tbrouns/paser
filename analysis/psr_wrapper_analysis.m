function psr_wrapper_analysis(filenames,loadPath,savePath,cfg)

% "cfg.fpath" should contain path to function

fpath = cfg.fpath;
names = split(fpath,'\');
names = names(~cellfun('isempty',names)); % remove empty cells
fname = names{end};
fpath = join(names(1:end-1),'\');

oldFolder = cd(fpath{1});
feval(fname,filenames,loadPath,savePath,cfg);
% try % TEMP COMMENT
%     feval(fname,filenames,loadPath,savePath,cfg);
% catch ME
%     disp(['Error in "' fname '"']);
%     disp(ME.message);
% end

cd(oldFolder);

end