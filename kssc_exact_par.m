function [ Z ] = kssc_exact_par( X, lambda, num_neighbors )

[D, N] = size(X);

if (num_neighbors >= N)
    num_neighbors = N - 1;
end

full_rows_mat = double(vl_kdtreequery(vl_kdtreebuild(X),X,X,'NUMNEIGHBORS',num_neighbors+1,'MAXNUMCOMPARISONS',10*(num_neighbors+1/10)));

[rows, cols] = process_neighbours(full_rows_mat);

inds_mat = full_rows_mat(2:end, :);

max_iterations = 200;

mu = 0.5;
mu_max = 10;

gamma_0 = 1.1;

tol_1 = 1*10^-3;
tol_2 = 1*10^-4;

Z_vals = zeros(num_neighbors, N);

parfor j = 1 : N
    
    Z_vals(:, j) = par_helper_exact(X(:, j), X(:, inds_mat(:,j)), lambda, D, num_neighbors, max_iterations, mu_max, mu, gamma_0, tol_1, tol_2);

end

Z = sparse(rows, cols, Z_vals(:), N, N);

end

