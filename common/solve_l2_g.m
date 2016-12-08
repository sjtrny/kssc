function [ x ] = solve_l2_g( w, lambda )
%MYSOLVE_L2 Summary of this function goes here
%   Detailed explanation goes here

% min lambda \lambda |x|_2 + 1/2 * |x-w|_2^2

% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end

end

