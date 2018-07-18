function psr_sst_ost_setpath(basepathCode)

%--below, no changes necessary
path(path,[basepathCode ]);
path(path,[basepathCode '/code/continuous/']);
path(path,[basepathCode '/code/continuous/neuralynx']);
path(path,[basepathCode '/code/continuous/txt']);
path(path,[basepathCode '/code/sortingNew/']);
%path(path,[basepathCode '/code/sortingNew/noiseRemoval']);
path(path,[basepathCode '/code/sortingNew/projectionTest']);
path(path,[basepathCode '/code/sortingNew/model']);
path(path,[basepathCode '/code/sortingNew/model/detection']);
%path(path,[basepathCode '/code/sortingNew/evaluation']);
path(path,[basepathCode '/code/osortTextUI']); % text user interface
path(path,[basepathCode '/code/osortGUI']); %graphical user interface
path(path,[basepathCode '/code/helpers']);
path(path,[basepathCode '/code/plotting']);
path(path,[basepathCode '/code/3rdParty/gabbiani']);

if ispc
    path(path,[basepathCode '/code/3rdParty/neuralynxWindows']);
else
    path(path,[basepathCode '/code/3rdParty/neuralynxUnixAll/binaries']);
    
end

path(path,[basepathCode '/code/3rdParty/cwtDetection']);

end