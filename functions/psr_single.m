function x = psr_single(x,parameters)

precision = 10^round(parameters.general.precision);
x = single(x) / precision;

end