function [ X ] = gen_unilrr_data( dim_data, dim_space, cluster_size, n_space, lower, upper )
%GEN_SUBSPACE_DATA Summary of this function goes here
%   Detailed explanation goes here


    X = zeros(dim_data, cluster_size*n_space);

    [T, ~] = qr(rand(dim_data,dim_data));
    if det(T) < 1
        T(:,1) = -T(:,1);
    end
    
    U = orth(rand(dim_data, dim_space));
    
    X(:, 1:cluster_size) = U * (lower + (abs(upper - lower)* rand(dim_space, cluster_size)));
    for j = 1 : n_space-1
        U = T * U;
        X(:, cluster_size*j+1: cluster_size*(j+1)) = U * (lower + (abs(upper - lower)* rand(dim_space, cluster_size)));
    end

end