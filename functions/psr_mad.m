function y = psr_mad(x,DIM)

% PSR_MAD - Calculate median absolute deviation (MAD)
%
% Syntax:  y = psr_mad(x,DIM)
%
% Inputs:
%    x   - Data array
%    DIM - Dimension along which we want to calculate the MAD
%
% Outputs:
%    y - The MAD value

% PASER: Processing and Analysis Schemes for Extracellular Recordings 
% https://github.com/tbrouns/paser

% Author: Terence Brouns
% Radboud University, Neurophysiology Dept. 
% E-mail address: t.s.n.brouns@gmail.com
% Date: 2018

%------------- BEGIN CODE --------------

if (nargin < 2)
    if (size(x,1) > size(x,2)); DIM = 1;
    else,                       DIM = 2;
    end
end

y = median(abs(x),DIM) / 0.6745;

end
