function [ Z ] = kssc_exact_binned_par( X, lambda, num_neighbors, n_threads )

[D, N] = size(X);

if (num_neighbors >= N)
    num_neighbors = N - 1;
end

full_rows_mat = double(vl_kdtreequery(vl_kdtreebuild(X),X,X,'NUMNEIGHBORS',num_neighbors+1,'MAXNUMCOMPARISONS',10*(num_neighbors+1/10)));

[rows, cols] = process_neighbours(full_rows_mat);

inds_mat = full_rows_mat(2:end, :);

max_iterations = 200;

mu = 0.1;
mu_max = 1;

gamma_0 = 1.1;

tol_1 = 1*10^-3;
tol_2 = 1*10^-4;



z_column_cells = cell(n_threads, 1);
ind_column_cells = cell(n_threads, 1);

inds = 1 : N;

n_cols = floor(N / n_threads);

for k = 1 : n_threads - 1
    ind_column_cells{k} = inds(:, (k-1)*n_cols + 1 : k*n_cols);
end
ind_column_cells{n_threads} = inds(:, (n_threads-1)*n_cols + 1 : N);


parfor h = 1 : n_threads
    
    thread_N = size( ind_column_cells{h}, 2);
    thread_Z = zeros(num_neighbors, thread_N);
    thread_inds = ind_column_cells{h};
    
    for j = 1 : thread_N
        cur_ind = thread_inds(j);
        
        thread_Z(:, j) = par_helper(X(:, cur_ind), X(:, inds_mat(:,cur_ind)), lambda, D, num_neighbors, max_iterations, mu_max, mu, gamma_0, tol_1, tol_2);
    end

    z_column_cells{h} = thread_Z;
    
end

Z_vals = zeros(num_neighbors, N);

for k = 1 : n_threads - 1
    Z_vals(:, (k-1)*n_cols + 1 : k*n_cols) = z_column_cells{k}; 
end
Z_vals(:, (n_threads-1)*n_cols + 1 : N) = z_column_cells{n_threads};

Z = sparse(rows, cols, Z_vals(:), N, N);

end

