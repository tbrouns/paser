function psr_remove_path(pathDel)

pathsAll = path;
pathsAll = strsplit(pathsAll,';');
k = ~cellfun(@isempty,strfind(pathsAll,pathDel)); %#ok
pathsDel = pathsAll(1,k);
pathsDel = strjoin(pathsDel,';');
rmpath(pathsDel);

end