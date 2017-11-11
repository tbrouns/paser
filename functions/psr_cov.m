% get covariance matrix of background nosie by randomly sampling 10000 timepoints
function c = psr_cov(data, samples)

nSamples = size(data,1);
nChans   = size(data,2);

max_samples = 10000;
waves       = zeros([max_samples samples nChans]);
data_index  = ceil((nSamples-samples) .* rand([1 max_samples]));

for j = 1:max_samples
    waves(j,:,:) = data(data_index(j) + (0:samples-1),:); % use bsxfun instead
end

c = cov(waves(:,:));