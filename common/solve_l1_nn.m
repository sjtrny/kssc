function [ x ] = solve_l1_nn( b, lambda )
% Function to solve soft thresholding problem
%
% arg min_{x} ||x - b||_{2}^{2} + lambda*||x||_{1}
%

x = max( sign(b).*max(abs(b) - lambda, 0), 0);

end

