% What do we want to test?
% How many true positives does knn choose
% How many false positives does knn choose
%
% This is dependant on ambient dimension, subspace dimension, distribution in the subspace 
% Also dependant on K in KNN

% 2D line plot (multiple with other parameters being the multiplew)
% X axis: increasing K
% Y axis: increasing true positives - increasing false positives

% 2D color grid plot
% X axis: increasing K
% Y axis: ambient/subspace dimension and distribution
% Grid color: true positives - false positives

% RESULTS:
% Strongest dependancy on mean and variance. In other words the structure
% within each subspace greatly affects success rate.
% Next greatest effect is the cluster size
% Then the ambient dimension


paths = ['common:'];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

% rng(1);

ambient_dim = 100;
sub_dim = 4;
n_space = 5;
cluster_size = 20;
% sub_mean = 0.5 ;
sub_mean = 0.05;
sub_var = 0.01;
% sub_var = 0.00025;
% sub_var = 0.5;

% ambient_dim = 32256;
% sub_dim = 9;
% n_space = 5;
% cluster_size = 20;
% sub_mean = 0.1;
% % sub_var = 0.00025;
% sub_var = 0.01;


n_samples = n_space * cluster_size;

cluster_inds = reshape(repmat(1:n_space, cluster_size, 1), 1, n_samples );

X = gen_depmultivar_data(ambient_dim, sub_dim, cluster_size, n_space, sub_mean, sub_var);

% X = normalize(X);

true_pos = get_block_diag(n_samples, n_space);
false_pos = 1 - true_pos;

max_nn = 100;

true_pos_rate = zeros(max_nn - 1,1);
false_pos_rate = zeros(max_nn - 1,1);

% figure, quiver3(0,0,0, 1, 0, 0)
% hold on, quiver3(0,0,0, 0, 1, 0)
% hold on, quiver3(0,0,0, 0, 0, 1)
% 
% line_specs = {'-ok', '-or', '-og', '-ob', '-oc'};
% line_col = {'k', 'r', 'g', 'b', 'c'};
% 
% hold on
% scatter3(X(1,1:20), X(2,1:20), X(3,1:20), line_col{1});
% 
% for k = 1 : 4
%     scatter3(X(1,20*k+1: 20*(k+1)), X(2,20*k+1: 20*(k+1)), X(3,20*k+1: 20*(k+1)), line_col{k+1});
% end


for k = 1 : max_nn - 1;
    
    full_rows_mat = double(vl_kdtreequery(vl_kdtreebuild(X),X,X,'NUMNEIGHBORS',k+1,'MAXNUMCOMPARISONS',10*(k+1/10)));
    rows_mat = reshape(full_rows_mat(2:end, :), k*n_samples, 1);


%     Z = zeros(n_samples, n_samples);
%     full_rows_mat = zeros(k, n_samples);
% 
%     % Build distance matrix
%     for i = 1 : n_samples
%         Z(i, :) = sum(abs(repmat(X(:,i), 1, n_samples) - X).^2).^(1/2);
%         [b, index] = sort(Z(i,:));
%         full_rows_mat(:, i) = index(1:k);
%     end
% 
%     rows_mat = reshape(full_rows_mat, k*n_samples, 1);
    
    cols_mat = reshape(repmat(1:n_samples, k, 1), k*n_samples, 1);
    
    index_mat = sparse(rows_mat, cols_mat, ones(k*n_samples, 1), n_samples, n_samples);
    
    true_pos_rate(k) = nnz(index_mat .* true_pos);
    false_pos_rate(k) = nnz(index_mat .* false_pos);
    
end

normalized_true_pos = true_pos_rate ./ ((1:(max_nn-1))'*n_samples);
normalized_false_pos = false_pos_rate ./ ((1:(max_nn-1))'*n_samples);

enough_true_pos = (true_pos_rate ./ (n_samples)) >= sub_dim;

figure, plot(normalized_true_pos)
% figure, plot(normalized_false_pos)
hold on
plot(enough_true_pos, 'g');

rmpath(paths);