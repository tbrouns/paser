function y = psr_sst_normpdf_stability(x,lambda,parameters)

y = normpdf(x,lambda,lambda^(parameters.cluster.stability.fvar));
y = y / sum(y);

end