paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end

% rng(1);

rows = 50;
n_space = 5;
cluster_size = 20;
subspace_dim = 10;

tic;
A = gen_depmultivar_data(rows, subspace_dim, cluster_size, n_space, 0.1, 2);
toc;

% tic;
% 
% A = zeros(rows, n_space *cluster_size);
% 
% for t = 1 : n_space
% 
%     basis = orth(randn(rows, subspace_dim));
% 
%     A(:, cluster_size*(t-1) + 1: cluster_size*t) = basis * ( rand(subspace_dim, cluster_size) - 0.5);
% 
% end
% 
% toc;

% max noise multiplier ~ 0.05

B = normalize(A);

w = randn(size(A)) * 0.05;
X = B + w;

% X = A;

X_normed = normalize(X);

truth = reshape(repmat(1:n_space, cluster_size, 1), n_space * cluster_size, 1);


% load('data/TIRlib');
% 
% n_spectra = 3;
% n_cluster = 5;
% cluster_size = 1000;
% 
% B = zeros(321, n_cluster * cluster_size);
% 
% for j = 1 : n_cluster
%     spectra_indices = randi(size(A,2),n_spectra,1);
% 
%     B(:,((j-1)*cluster_size)+1:j*cluster_size) = A(:, spectra_indices) * gen_depcoeff(0.1, 0.001, 3, cluster_size);
% end
% 
% truth = reshape(repmat(1:n_cluster, cluster_size, 1), n_cluster * cluster_size, 1);
% 
% w = randn(size(B)) * 0.05;
% X = B + w;
% 
% X_normed = normalize(B);

% SSC
% 
% tic;
% [Z_ssc_exact] = ssc_exact(X_normed, 0.1);
% toc;

% tic;
% [Z_ssc_exact_par] = ssc_exact_par(X_normed, 0.1);
% toc;
% 
% tic;
% [Z_ssc_relaxed] = ssc_relaxed(X_normed, 0.1);
% toc;
% 
tic; 
[Z_ssc_relaxed_par] = ssc_relaxed_par(X_normed, 0.1);
toc;

[ssc_clusters,~,~] = ncutW((abs(Z_ssc_relaxed_par)+abs(Z_ssc_relaxed_par')), n_space);
clusters_ssc = condense_clusters(ssc_clusters,1);
missrate_ssc = Misclassification(clusters_ssc, truth);



% Scalable SSC

% tic;
% [Z_kssc_exact] = kssc_exact(X_normed, 0.1, 25);
% toc;
% 
% tic;
% [Z_kssc_exact_par] = kssc_exact_par(X_normed, 0.1, 25);
% toc;
% 
% tic;
% [Z_kssc_relaxed] = kssc_relaxed(X_normed, 0.1, 25);
% toc;

% tic;
% [Z_kssc_relaxed_par] = kssc_relaxed_par(X_normed, 0.1, 25);
% toc;

% [kssc_clusters,~,~] = ncutW((abs(Z_kssc_relaxed_par)+abs(Z_kssc_relaxed_par')), n_space);
% clusters_kssc = condense_clusters(kssc_clusters,1);
% missrate_kssc = Misclassification(clusters_kssc, truth);

% tic;
% Z_tsc = tsc_sc(X_normed, 20);
% toc;
% 
% tic;
% clusters_sssc = get_sssc_clusters(X_normed, 20, 0.1, n_space );
% toc;
% 
% tic;
% Z_gsc = gsc_sc(X_normed, 10, 5);
% toc;
% 
% [gsc_clusters,~,~] = ncutW((abs(Z_gsc)+abs(Z_gsc')), n_space);
% clusters_gsc = condense_clusters(gsc_clusters,1);
% missrate_gsc = Misclassification(clusters_gsc, truth);
