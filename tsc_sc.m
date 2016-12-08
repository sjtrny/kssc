function [ Z ] = tsc_sc( X, num_neighbours )
%TSC_SC Summary of this function goes here
%   Detailed explanation goes here

[~, n] = size(X);

full_rows_mat = double(vl_kdtreequery(vl_kdtreebuild(X),X,X,'NUMNEIGHBORS',num_neighbours+1,'MAXNUMCOMPARISONS',10*(num_neighbours+1/10)));

[inds_rows, inds_cols] = process_neighbours(full_rows_mat);

inds_rows = reshape(inds_rows, num_neighbours, n);

Z_vals = zeros(num_neighbours, n);

parfor k = 1 : n
    
    in_prod = X(:, inds_rows(:, k))' *  X(:,k);
    in_prod(in_prod > 1) = 1;
    
    Z_vals(:,k) = exp( -2 * acos(in_prod) );
    
end

Z = sparse(inds_rows(:), inds_cols(:), Z_vals(:), n, n);


% [~, n] = size(X);
% 
% Z_vals = zeros(num_neighbours, n);
% rows = zeros(num_neighbours, n);
% 
% parfor k = 1 : n
%     sub_X = X;
%     sub_X(:,k) = [];
%     
%     inds = double(vl_kdtreequery(vl_kdtreebuild(sub_X),sub_X, X(:,k) ,'NUMNEIGHBORS',num_neighbours,'MAXNUMCOMPARISONS',10*(num_neighbours/10)));
%     
%     in_prod = sub_X(:, inds)' *  X(:,k);
%     in_prod(in_prod > 1) = 1;
%     
%     Z_vals(:,k) = exp( -2 * acos(in_prod) );
%     rows(:,k) = inds;
%     
% end
% 
% Z = sparse(rows(:), reshape( repmat(1:n, num_neighbours, 1) , num_neighbours * n, 1), Z_vals(:), n, n);

end

