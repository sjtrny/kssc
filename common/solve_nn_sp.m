function [ X, s ] = solve_nn_sp( Y, tau, rows, cols )
%NN_PROX
%   This function solves the proximal nuclear norm problem
% 
%   min lambda * |X|_* + 1/2*|X - Y|^2
%
%   solved by singular value thresholding
% 
%   Written by Stephen Tierney

[U, S, V] = svds(Y);

s = diag(S);

ind = find(s <= tau);
s(ind) = 0;

ind = find(s > tau);
s(ind) = s(ind) - tau;

S = diag(s);

if (size(Y,1) ~= size(Y,2))
    rows = size(Y,1);
    cols = size(Y,2);
    S(:,rows+1:cols) = zeros(rows,cols - rows);
end

n = size(Y, 1);

X = sparse(rows, cols, sp_mult(U*S, V', rows, cols), n, n);

% X = U*S*(V');

% imagesc(X-X_new), colorbar;

end

