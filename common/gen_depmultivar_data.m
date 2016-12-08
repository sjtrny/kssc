function [ X ] = gen_depmultivar_data( dim_data, dim_space, cluster_size, n_space, m, v )
%% Generate Dependant Multivariate Data
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

    X = zeros(dim_data, cluster_size*n_space);

    for j = 1 : n_space
        
        basis = orth(randn(dim_data, dim_space));
        
%         X(:, cluster_size*(j-1) + 1: cluster_size*j) = basis * gen_depcoeff(m, v, dim_space, cluster_size);

        X(:, cluster_size*(j-1) + 1: cluster_size*j) = basis * (mvnrnd( m * ones(1, dim_space), v * ones(1,dim_space, cluster_size))');
        
    end


end

