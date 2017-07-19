files = dir('Spikes*.mat');
files = char(files.name);
numfiles = length(files(:,1));

for iFile = 1:numfiles
    load(files(iFile,:));
    ums_visualization(spikes,metadata)
end

