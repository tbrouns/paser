function x = psr_int16_to_single(x,parameters)
if (isa(x,'int16'))
    precision = 10^parameters.general.precision;
    x = single(x) / precision;
end
end