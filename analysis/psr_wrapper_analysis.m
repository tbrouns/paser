function psr_wrapper_analysis(cfg)

% "cfg.fpath" should contain path to function

fpath = cfg.fpath;
names = split(fpath,'\');
names = names(~cellfun('isempty',names)); % remove empty cells
fname = names{end};
fpath = join(names(1:end-1),'\');

oldFolder = cd(fpath{1}); % cd to path to ensure we call correct function
feval(fname,cfg);
cd(oldFolder); % go back to where we came from

end