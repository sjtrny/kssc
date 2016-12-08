function [ Z ] = ssc_exact_par( X, lambda )

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

parfor j = 1 : N
    
    func_vals = zeros(max_iterations, 1);
    
    x_cur = X(:, j);
    
    z = zeros(N, 1);
    
    e = zeros(D, 1);
    
    y = zeros(D, 1);
    
    mu = 0.1;
    
    for k = 1 : max_iterations

        z_prev = z;
        e_prev = e;

        % Update Z
        z = solve_l1_nn(z - 1/rho * (mu*X'*(X*z - x_cur + e + 1/mu * y)), lambda/rho);
        z(j) = 0;

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
    
    Z(:, j) = z;
    
end



end

