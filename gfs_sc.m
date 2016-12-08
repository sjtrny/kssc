function [ Z ] = gfs_sc( X, num_neighbors )
%GFS_SC Summary of this function goes here
%   Detailed explanation goes here
% 

[m, n] = size(X);

id_mat = eye(m, m);

Z_vals = zeros(num_neighbors, n);
i_inds_out = zeros(num_neighbors, n);

global_inds = 1:n;

for k = 1 : n
    
    i_inds = zeros(num_neighbors, 1);
    
    local_X = zeros(m, num_neighbors);
    
    y = X(:, k);
    
    s = y;
    
    sub_X = X;
    sub_X(:,k) = [];

    global_inds_map = global_inds;
    global_inds_map(:,k) = [];
    
    normR = Inf;
    
    j = 0;
    
    while (j < num_neighbors && normR > 1*10^-3)
        
        j = j + 1;
        
%         kdtree = vl_kdtreebuild(sub_X);
%         local_ind = double(vl_kdtreequery(kdtree, sub_X, s, 'NUMNEIGHBORS', 1, 'MAXNUMCOMPARISONS', 10*(1/10)));
        
        [~, local_ind] = max(sub_X' * s);

        i_inds(j) = global_inds_map(local_ind);
        
        local_X(:, j) = X(:, i_inds(j));
        
        s = (id_mat - local_X(:,1:j) * pinv(local_X(:,1:j))) * y;

%         x_j = local_X(:,1:j)' * local_X(:,1:j) \ local_X(:,1:j)' * y;
%         s = y - (local_X(:,1:j) * x_j);
        
%         sub_X = X;
%         sub_X(:,[k, i_inds(j)]) = [];
%         
%         global_inds_map = global_inds;
%         global_inds_map(:,[k, i_inds(j)]) = [];
%         
        normR = norm(s, 2);
        
    end

%     Z_vals(:, 1:j) = 1;
    
%     Z_vals(:,  k) = pinv(local_X) * y; 

%     Z_vals(1:j, k) = local_X(:,1:j)' * local_X(:,1:j) \ local_X(:,1:j)' * y;
%     Z_vals(:, k) = local_X(:,1:j)' * local_X(:,1:j) \ local_X(:,1:j)' * y;
    Z_vals(:, k) = pinv(local_X(:,1:j)) * y;

    i_inds_out(:, k) = i_inds;
    
end

cols = reshape(repmat(1:n, num_neighbors, 1), n*num_neighbors, 1);

inds_vec = i_inds_out(:);

[a, ~] = find(inds_vec);

Z = sparse(inds_vec(a), cols(a), ones(length(a), 1), n, n);


end

