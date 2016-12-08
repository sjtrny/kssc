function [ Z ] = kssc_relaxed_par( X, lambda, num_neighbors )
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda || Z ||_1
%
% by accelerated gradient descent
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 200;

[D, N] = size(X);

if (num_neighbors >= N)
    num_neighbors = N - 1;
end

full_rows_mat = double(vl_kdtreequery(vl_kdtreebuild(X),X,X,'NUMNEIGHBORS',num_neighbors+1,'MAXNUMCOMPARISONS',10*(num_neighbors+1/10)));

[rows, cols] = process_neighbours(full_rows_mat);

inds_mat = full_rows_mat(2:end, :);

Z_vals = zeros(num_neighbors, N);

rho = 1;
gamma = 1.1;

alpha = 1;

tol = 1*10^-6;

parfor j = 1 : N
    
    Z_vals(:, j) = par_helper_relaxed(X(:, j), X(:, inds_mat(:,j)), lambda, D, num_neighbors, max_iterations, rho, gamma, alpha, tol);
    
end

Z = sparse(rows, cols, Z_vals(:), N, N);



end

