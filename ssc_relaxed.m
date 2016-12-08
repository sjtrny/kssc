function [ Z ] = ssc_relaxed( X, lambda )
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

N = size(X, 2);

Z = zeros(N);
J = zeros(N);

rho = 1;
gamma = 1.1;

alpha = 1;

diag_inds = ((1:N) * N) + (1:N) - N;

for k = 1 : max_iterations

    Z_prev = Z;
    alpha_prev = alpha;
    
    searching = true;
    while( searching )
        partial = -X'*(X - X*J);

        Z = solve_l1_nn(J - 1/rho * partial, lambda/rho);
        
        Z(diag_inds) = 0;
    
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
    
    
    if ( abs(func_vals(k) - previous_func_val) <= 1*10^-6 )
        break;
    else
        previous_func_val = func_vals(k);
    end

end

end

