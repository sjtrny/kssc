function [ Z ] = kssc_relaxed( X, lambda, num_neighbors )
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

func_vals = zeros(max_iterations, 1);
previous_func_val = Inf;

[D, N] = size(X);

if (num_neighbors >= N)
    num_neighbors = N - 1;
end

full_rows_mat = double(vl_kdtreequery(vl_kdtreebuild(X),X,X,'NUMNEIGHBORS',num_neighbors+1,'MAXNUMCOMPARISONS',10*(num_neighbors+1/10)));

[rows, cols] = process_neighbours(full_rows_mat);

% inds_mat = full_rows_mat(2:end, :);

Z = sparse(N, N);
J = sparse(N, N);

rho = 1;
gamma = 1.1;

alpha = 1;

tol = 1*10^-6;


for k = 1 : max_iterations

    Z_prev = Z;
    alpha_prev = alpha;
    
    searching = true;
    while( searching )
        partial = sparse(rows, cols, selmult_mt_m_fast(-X, X - X*J,  rows, cols), N, N);
        
        Z = solve_l1_spm(J - 1/rho * partial, lambda/rho);
        
%         partial2 = -X' * (X - X*full(J));
%         Z2 = solve_l1_nn(full(J) - 1/rho * full(partial2), lambda/rho);
        
        func_vals(k, 1) = lambda*norm_l1(Z) + 1/2 *norm(X - X*Z, 'fro')^2;
        
        approx = lambda*norm_l1(Z) + 0.5*norm(X - X*J, 'fro')^2 + sum(sum((Z - J).*partial)) + 0.5*rho*norm(Z - J, 'fro')^2;

        if ( func_vals(k, 1) > approx )
            rho = gamma * rho;
        else
            searching = false;
        end
    end
    
    alpha = (1 + sqrt(1 + 4*alpha_prev^2)) / 2;
    
    J = Z + ((alpha_prev - 1)/alpha)*(Z - Z_prev);
    
    
    if ( abs(func_vals(k) - previous_func_val) <= tol )
        break;
    else
        previous_func_val = func_vals(k);
    end

end

[r, c] = ind2sub([N, N], find(Z));
Z = sparse(r, c, Z(find(Z)), N, N);

end

