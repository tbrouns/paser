function data = psr_lfp_mean_subtraction(data)

nBlocks = length(data);
for iBlock = 1:nBlocks
    data{iBlock}.trial{1} = meanSubtraction(data{iBlock}.trial{1});
end

end

function data = meanSubtraction(data)

dataMean = mean(data,1);
data = bsxfun(@minus,data,dataMean);

end