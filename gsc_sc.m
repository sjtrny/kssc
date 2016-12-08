function [ Z ] = gsc_sc( X, num_neighbors, max_sub_dim )
%GSC_SC Summary of this function goes here
%   Detailed explanation goes here

[~, n] = size(X);

rows = zeros(num_neighbors + 1, n);

for i = 1 : n
    
    inds_set = zeros(num_neighbors + 1, 1);
    inds_set(1) = i;
    
    U = [];
    
    for k = 1 : num_neighbors
        
        if k < max_sub_dim
            U = X(:, inds_set(1:k));
        end
        
        remaining_inds = 1:n;
        
        remaining_inds(inds_set(1:k)) = [];
        
%         remaining_inds = setdiff( 1:n , inds_set(1:k));
        
        in_prods = U' * X(:, remaining_inds);
        
        [~, max_ind] = max( sqrt( sum(in_prods .^ 2, 1) ) );

        inds_set(1+k) = remaining_inds(max_ind);
        
    end
    
    rows(:, i) = inds_set(1:end);
    
end

cols = reshape( repmat(1:n, num_neighbors + 1, 1) , (num_neighbors+1) * n, 1);

Z = sparse(rows(:) , cols , ones((num_neighbors+1) * n, 1), n, n);

end

