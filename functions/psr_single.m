function x = psr_single(x,parameters)
precision = 10^parameters.general.precision;
x = single(x) / precision;
end