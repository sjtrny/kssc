function [ Z ] = ssc_exact_binned_par( X, lambda, n_threads )

max_iterations = 200;

N = size(X, 2);
D = size(X, 1);

mu_max = 1;
rho = (norm(X,2)^2) * 1.2;
gamma_0 = 1.1;

tol_1 = 1*10^-3;
tol_2 = 1*10^-4;

normfX = norm(X,'fro');

Z = zeros(N, N);

z_column_cells = cell(n_threads, 1);
x_column_cells = cell(n_threads, 1);
ind_column_cells = cell(n_threads, 1);

inds = 1 : N;

n_cols = floor(N / n_threads);

for k = 1 : n_threads - 1

    x_column_cells{k} = X(:, (k-1)*n_cols + 1 : k*n_cols);
    ind_column_cells{k} = inds(:, (k-1)*n_cols + 1 : k*n_cols);
    
end

x_column_cells{n_threads} = X(:, (n_threads-1)*n_cols + 1 : N);
ind_column_cells{n_threads} = inds(:, (n_threads-1)*n_cols + 1 : N);

parfor h = 1 : n_threads
    
    thread_N = size( x_column_cells{h}, 2);
%     thread_X = x_column_cells{h};
    thread_Z = zeros(N, thread_N);
    thread_inds = ind_column_cells{h};
    
    for j = 1 : thread_N
        
        func_vals = zeros(max_iterations, 1);
    
%         x_cur = thread_X(:, j);
        cur_ind = thread_inds(j);
        x_cur = X(cur_ind);

        z = zeros(N, 1);

        e = zeros(D, 1);

        y = zeros(D, 1);

        mu = 0.1;
    
        for k = 1 : max_iterations

            z_prev = z;
            e_prev = e;

            % Update Z
            z = solve_l1_nn(z - 1/rho * (mu*X'*(X*z - x_cur + e + 1/mu * y)), lambda/rho);
            
            z(cur_ind) = 0;

            % Update E
            e = solve_l2(x_cur - X*z - 1/mu * y, 1/mu);

            % Check convergence
            func_vals(k) = 0.5*norm(e, 2)^2 + lambda*norm_l1(z);

            stop_2 = mu * max([ sqrt(rho)*norm(z - z_prev,2), norm(e - e_prev, 2)]) / normfX;

            if ( (norm(X*z - x_cur + e, 'fro')/normfX) < tol_1 ...
                && stop_2 < tol_2)
                break;
            end

            % Update Y
            y = y + mu*(X*z - x_cur + e);

            % Uppdate rate parameters
            if (stop_2 < tol_2)
                gamma = gamma_0;
            else
                gamma = 1;
            end

            mu = min(mu_max, gamma * mu);

        end
        
        thread_Z(:, j) = z;
        
        
    end
    
    z_column_cells{h} = thread_Z;
    
end

for k = 1 : n_threads - 1

    Z(:, (k-1)*n_cols + 1 : k*n_cols) = z_column_cells{k};
end

Z(:, (n_threads-1)*n_cols + 1 : N) = z_column_cells{n_threads}; 


end