function data = psr_lfp_artifact_channel(data,parameters)

dataAbs = psr_lfp_conversion(data);
dataAbs = cat(2,dataAbs{:}); % Combine data across all blocks

nChans  = size(dataAbs,1);
dataAbs =  abs(dataAbs);
dataAvg = mean(dataAbs,2); 
dataStd =  std(dataAbs,[],2);

chanDiff = zeros(nChans,1);

for iChan = 1:nChans
   
    I = true(nChans,1);
    I(iChan) = false;
    
    dataStdChan = mean(dataStd(I));
    dataAvgMean = mean(dataAvg(I));
    
    chanDiff(iChan) = (dataAvg(iChan) - dataAvgMean) / dataStdChan;
    
end

% Remove channels
keep = chanDiff <= parameters.lfp.artifact.chan.thresh;
nBlocks = length(data);
for iBlock = 1:nBlocks
    data{iBlock}.trial{1} = data{iBlock}.trial{1}(keep,:);
    data{iBlock}.label    = data{iBlock}.label(keep);
end