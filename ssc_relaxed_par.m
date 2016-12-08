function [ Z ] = ssc_relaxed_par( X, lambda )
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

[~, N] = size(X);

Z = zeros(N, N);

gamma = 1.1;

tol = 1*10^-6;

parfor t = 1 : N
    
    func_vals = zeros(max_iterations, 1);
    previous_func_val = Inf;
    
    x_cur = X(:, t);
    z = zeros(N, 1);
    j = zeros(N, 1);
    rho = 1;
    alpha = 1;
    
    looping = true;
    k = 1;
    while(looping)
        
        z_prev = z;
        alpha_prev = alpha;

        searching = true;
        while( searching )
            partial = -X'*(x_cur - X*j);

            z = solve_l1_nn(j - 1/rho * partial, lambda/rho);

            z(t) = 0;

            func_vals(k, 1) = lambda*norm_l1(z) + 1/2 *norm(x_cur - X*z, 'fro')^2;

            approx = lambda*norm_l1(z) + 0.5*norm(x_cur - X*j, 'fro')^2 + sum((z - j).*partial) + 0.5*rho*norm(z - j, 'fro')^2;

            if ( func_vals(k, 1) > approx )
                rho = gamma * rho;
            else
                searching = false;
            end
        end
    
        alpha = (1 + sqrt(1 + 4*alpha_prev^2)) / 2;

        j = z + ((alpha_prev - 1)/alpha)*(z - z_prev);

        if ( abs(func_vals(k) - previous_func_val) <= tol ...
                || k >= max_iterations)
            looping = false;
        else
            previous_func_val = func_vals(k);
            k = k + 1;
        end
    end
    
    Z(:, t) = z;
    
end

end
