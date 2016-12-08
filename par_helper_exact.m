function [ z ] = par_helper_exact( x, X_mat, lambda, D, num_neighbors, max_iterations, mu_max, mu, gamma_0, tol_1, tol_2)
%ITER_HELPER Summary of this function goes here
%   Detailed explanation goes here

func_vals = zeros(max_iterations, 1);

rho = (norm(X_mat,2)^2) * 1.2;
normfx = norm(X_mat, 2);

z = zeros(num_neighbors, 1);
e = zeros(D, 1);
y = zeros(D, 1);


for k = 1 : max_iterations
        
    z_prev = z;
    e_prev = e;

    % Update z
    z = solve_l1(z_prev - 1/rho * (mu*X_mat' * (X_mat*z_prev - x + e + 1/mu * y)),  lambda/rho);
    
    % Update e
    e = solve_l2(x - X_mat*z - 1/mu * y, 1/mu);

    % Check convergence
    func_vals(k) = 0.5*norm(e,2)^2 + lambda*norm_l1(z);

    stop_2 = mu * max([ sqrt(rho)*norm(z - z,2), norm(e - e_prev, 2)]) / normfx;

    if ( (norm(X_mat*z - x + e, 'fro')/normfx) < tol_1 ...
        && stop_2 < tol_2)
        break;
    end

    % Update y
    y = y + mu*(X_mat*z - x + e);

    % Update rate parameters
    if (stop_2 < tol_2)
        gamma = gamma_0;
    else
        gamma = 1;
    end

    mu = min(mu_max, gamma * mu);

end


end

