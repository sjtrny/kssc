function [ z ] = par_helper_relaxed(x, X_mat, lambda, D, num_neighbors, max_iterations, rho, gamma, alpha, tol)
%ITER_HELPER Summary of this function goes here
%   Detailed explanation goes here

func_vals = zeros(max_iterations, 1);
previous_func_val = Inf;

z = zeros(num_neighbors, 1);
j = zeros(num_neighbors, 1);

for k = 1 : max_iterations

    z_prev = z;
    alpha_prev = alpha;
    
    searching = true;
    while( searching )
        partial = -X_mat'*(x - X_mat*j);

        z = solve_l1_nn(j - 1/rho * partial, lambda/rho);
    
        func_vals(k, 1) = lambda*norm_l1(z) + 1/2 *norm(x - X_mat*z, 'fro')^2;
        
        approx = lambda*norm_l1(z) + 0.5*norm(x - X_mat*j, 'fro')^2 + sum((z - j).*partial) + 0.5*rho*norm(z - j, 2)^2;

        if ( func_vals(k, 1) > approx )
            rho = gamma * rho;
        else
            searching = false;
        end
    end
    
    alpha = (1 + sqrt(1 + 4*alpha_prev^2)) / 2;
    
    j = z + ((alpha_prev - 1)/alpha)*(z - z_prev);

    if ( abs(func_vals(k) - previous_func_val) <= tol )
        break;
    else
        previous_func_val = func_vals(k);
    end

end


end

