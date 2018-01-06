function [B,MSEs] = jisotonic5(A,direction,weights)
% jisotonic5 - isotonic regression (jfm, may 2105)
%
% [B,MSEs] = jisotonic5(A,direction,weights)
%   A is the input a vector
%   direction = 'increading', 'decreasing', 'updown', or 'downup'
%   weights is the optional input weight vector (same size as A)
%   B is the output vector (same size as A)
%   MSEs is used internally for 'updown' and 'downup' directions
%
% Magland 5/19/2015

if (~isrow(A)); A = A'; end

if (nargin < 2); direction = 'increasing'; end
if (nargin < 3); weights = ones(size(A));  end

try
    [~,~] = jisotonic5_mex(1,1);
catch
    compile_mex_jisotonic5;
end

if (strcmp(direction,'decreasing'))
    [B,MSEs]=jisotonic5_mex(-A,weights); B=-B;
    return;
elseif (strcmp(direction,'updown'))
    [~,MSE1] = jisotonic5_mex(A,weights);
    [~,MSE2] = jisotonic5_mex(A(end:-1:1),weights(end:-1:1));
    MSE2 = MSE2(end:-1:1);
    MSE0 = MSE1 + MSE2;
    
    [~,best_ind] = min(MSE0);
    [C1,~] = jisotonic5_mex(A(1:best_ind),weights(1:best_ind));
    [C2,~] = jisotonic5_mex(-A(best_ind:end),weights(best_ind:end));
    C2 = -C2;
    B  = [C1(1:best_ind),C2(2:end)];
    if (isnan(B(1))); disp('ISO-SPLIT: jisotonic5: NaN'); end
    return;
elseif (strcmp(direction,'downup'))
    B = -jisotonic5(-A,'updown',weights);
    return;
else
    if (~strcmp(direction,'increasing'))
        disp(['ISO-SPLIT: invalid direction in jisotonic5: ' direction]);
        return;
    end
end

try
    [B,MSEs] = jisotonic5_mex(A,weights);
catch
    disp('Unable to run mex file -- You must compile using: mex jisotonic5_mex.cpp');
end

end