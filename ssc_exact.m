function [ Z, func_vals, k ] = ssc_exact( X, lambda )

max_iterations = 200;

func_vals = zeros(max_iterations, 1);

N = size(X, 2);

Z = zeros(size(X, 2));
E = zeros(size(X));
Y = zeros(size(X));

mu = 0.5;
mu_max = 1;
rho = (norm(X,2)^2) * 1.2;
gamma_0 = 1.1;

tol_1 = 1*10^-3;
tol_2 = 1*10^-4;

normfX = norm(X,'fro');
diag_inds = ((1:N) * N) + (1:N) - N;

for k = 1 : max_iterations

    Z_prev = Z;
    E_prev = E;
    
    % Update Z
    Z = solve_l1_nn(Z_prev - 1/rho * (mu*X'*(X*Z_prev - X + E + 1/mu * Y)), lambda/rho);
    Z(diag_inds) = 0;
    
    % Update E
    E = solve_l2(X - X*Z - 1/mu * Y, 1/mu);
    
    % Check convergence
    func_vals(k) = 0.5*norm(E, 'fro')^2 + lambda*norm_l1(Z);
    
    stop_2 = mu * max([ sqrt(rho)*norm(Z - Z_prev,'fro'), norm(E - E_prev, 'fro')]) / normfX;
    
    if ( (norm(X*Z - X + E, 'fro')/normfX) < tol_1 ...
        && stop_2 < tol_2)
        break;
    end
    
    % Update Y
    Y = Y + mu*(X*Z -X + E);
    
    % Uppdate rate parameters
    if (stop_2 < tol_2)
        gamma = gamma_0;
    else
        gamma = 1;
    end
    
    mu = min(mu_max, gamma * mu);
    
end

end

