function y = psr_mad(x)
    if (size(x,2) > size(x,1)); x = x'; end
    y = median(abs(x)) / 0.6745;
end
