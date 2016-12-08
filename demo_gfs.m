paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end

% rng(1);

rows = 100;
n_space = 5;
cluster_size = 20;
subspace_dim = 4;

tic;
A = gen_depmultivar_data(rows, subspace_dim, cluster_size, n_space, 0.1, 0.01);
toc;

B = normalize(A);

w = randn(size(A)) * 0.001;
X = B + w;

% X = B;

X_normed = normalize(X);

truth = reshape(repmat(1:n_space, cluster_size, 1), n_space * cluster_size, 1);

tic; 
[Z_gfs] = gfs_sc(X_normed, 10);
toc;

[gfs_clusters,~,~] = ncutW((abs(Z_gfs)+abs(Z_gfs')), n_space);
clusters_gfs = condense_clusters(gfs_clusters,1);
missrate_gfs = Misclassification(clusters_gfs, truth);


tic; 
[Z_ssc_relaxed_par] = ssc_relaxed_par(X_normed, 0.1);
toc;

[ssc_clusters,~,~] = ncutW((abs(Z_ssc_relaxed_par)+abs(Z_ssc_relaxed_par')), n_space);
clusters_ssc = condense_clusters(ssc_clusters,1);
missrate_ssc = Misclassification(clusters_ssc, truth);