function data = psr_lfp_mean_subtraction(data)

% Subtract channel mean from all probe channels

nBlocks = length(data);
for iBlock = 1:nBlocks
    dataBlock = data{iBlock};
    if (isfield(dataBlock,'trial'))
        dataBlock = dataBlock.trial{1};
        nChans = size(dataBlock,1);
        if (nChans > 1); data{iBlock}.trial{1} = meanSubtraction(dataBlock); end
    end
end

end

function data = meanSubtraction(data)

dataMean = nanmean(data,1); % Mean across channels
data = bsxfun(@minus,data,dataMean);

end