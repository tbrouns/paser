% get covariance matrix of background nosie by randomly sampling 10000 timepoints
function c = get_covs( data, samples )

    num_trials = length(data);
    num_channels = size(data{1},2);
    for j = 1:num_trials, num_samples(j) = size(data{j},1); end

    max_samples = 10000;
    waves = zeros( [max_samples samples num_channels] );
    tr_index = ceil( num_trials * rand([1 max_samples]) );
    data_index = ceil( (num_samples(tr_index)-samples) .* rand([1 max_samples]) );
    for j = 1:max_samples
       waves(j,:,:) = data{tr_index(j)}(data_index(j)+[0:samples-1],:);  
    end
        
   c = cov( waves(:,:) );
   
end