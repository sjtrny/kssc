function [ Z, func_vals, k ] = ssc_exact_old( X, lambda )

max_iterations = 200;

func_vals = zeros(max_iterations, 1);

Z = zeros(size(X, 2));

E = zeros(size(X));

Y = zeros(size(X));

mu = 0.1;
mu_max = 1;
rho = (norm(X,2)^2) * 1.2;

gamma_0 = 1.1;

covar = X'*X;

normfX = norm(X,'fro');

tol_1 = 1*10^-3;
tol_2 = 1*10^-4;

diag_inds = find(speye(size(Z)));

for k = 1 : max_iterations

    Z_prev = Z;
    E_prev = E;
    
    % Solve for Z
    
    partial = mu*(covar*Z - covar - X'*(- E - 1/mu * Y));
    
    V = Z - 1/rho * partial;
    
    Z = solve_l1(V, lambda/rho);
    
    Z(diag_inds) = 0;
    
    Z(Z < 0) = 0;
    
    % Solve for E
    
    V = X - X*Z - 1/mu * Y;
    
    E = solve_l2(V, 1/mu);
    
    % Update Y
    
    Y = Y + mu*(X*Z -X + E);
    
    if (mu * max(sqrt(rho) * norm(Z - Z_prev,'fro'), norm(E - E_prev))/norm(X,'fro') < tol_2)
        gamma = gamma_0;
    else
        gamma = 1;
    end
    
    mu = min(mu_max, gamma * mu);
    
    % Check convergence
    
    func_vals(k) = 0.5*norm(E, 'fro')^2 + lambda*norm_l1(Z);
    
    if ( (norm(X*Z - X + E, 'fro')/normfX < tol_1) ...
        && (mu * max([ sqrt(rho)*norm(Z - Z_prev,'fro'), norm(E - E_prev, 'fro')]) / normfX < tol_2))
        break;
    end
    
%     if ( abs(func_vals(k) - previous_func_val) <= 1*10^-6 )
%         break;
%     else
%         previous_func_val = func_vals(k);
%     end
    
end

end

