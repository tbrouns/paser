function y = psr_mad(x,DIM)

    if (nargin < 2)
        if (size(x,1) > size(x,2)); DIM = 1; 
        else,                       DIM = 2;
        end
    end
    
    y = median(abs(x),DIM) / 0.6745;
    
end
