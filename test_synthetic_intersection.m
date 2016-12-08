paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

rng(1);

rows = 200;

subspace_dim = 10;

n_cores = 8;

if ~exist('matlabpool','file') && ~matlabpool('size') == 0
    matlabpool close
end

matlabpool(n_cores);

n_steps = 10;
n_runs = 20;

ssc_accuracy = zeros(n_steps + 1, n_runs);
kssc_accuracy = zeros(n_steps + 1, n_runs);

ssc_runtime = zeros(n_steps + 1, n_runs);
kssc_runtime = zeros(n_steps + 1, n_runs);

m = 0.01;
v = 1;

n_space = 2;

for k = 1 : n_steps + 1
    
    cluster_size = 20 * subspace_dim;
    
    truth = reshape(repmat(1:n_space, cluster_size, 1), n_space *cluster_size, 1);
    
    for i = 1 : n_runs
            
        s1_basis = orth(randn(rows, subspace_dim));

        % s - number of intersecting dimensions
        s = k - 1;

        s2_basis_start = orth(randn(rows, n_steps - s));

        s2_basis = [s2_basis_start s1_basis(:, subspace_dim - s + 1 : subspace_dim)];

        A(:, 1 : cluster_size) = s1_basis * rand(subspace_dim, cluster_size);
        A(:, cluster_size+1 : cluster_size*2) = s2_basis * rand(subspace_dim, cluster_size);

        X_normed = normalize(A);

        tic;
        [Z_ssc_relaxed] = ssc_relaxed(X_normed, 0.1);
        ssc_runtime(k, i) = toc;

        clusters_ssc = condense_clusters(ncutW((abs(Z_ssc_relaxed)+abs(Z_ssc_relaxed')), n_space),1);
        ssc_accuracy(k, i) = Misclassification(clusters_ssc, truth);

        tic;
        [Z_kssc_relaxed] = kssc_relaxed_par(X_normed, 0.1, 10);
        kssc_runtime(k, i) = toc;

        clusters_kssc = condense_clusters(ncutW((abs(Z_kssc_relaxed)+abs(Z_kssc_relaxed')), n_space),1);
        kssc_accuracy(k, i) = Misclassification(clusters_kssc, truth);
    
    end
    
end

matlabpool close


h1 = plot(mean(ssc_accuracy,2), '-*m');
hold on
h2 = plot(mean(kssc_accuracy,2), '-*k');

set(gca, 'fontsize', 14);

xlim([1, 10]);

legend([h1, h2], 'SSC', 'kSSC', 'Location','northwest')

xlabel('Dimension of Intersection', 'FontSize', 18);

ylabel('Subspace Clustering Error (SCE)', 'FontSize', 18);

print(gcf, '-depsc2', 'test_synthetic_intersection.eps');

close all