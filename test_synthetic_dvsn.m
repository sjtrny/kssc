paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

rng(1);

n_cores = 8;

if ~exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool close
end

matlabpool(n_cores);

n_steps = 10;
n_runs = 50;

ambient_list = [30, 50, 100];

ssc_accuracy = zeros(n_steps, n_steps, n_runs, length(ambient_list));
kssc_accuracy = zeros(n_steps, n_steps, n_runs, length(ambient_list));

ssc_runtime = zeros(n_steps, n_steps, n_runs, length(ambient_list));
kssc_runtime = zeros(n_steps, n_steps, n_runs, length(ambient_list));

n_space = 5;

max_cluster_size = 150;
max_dimension = 30;

for k = 1 : n_steps

    cluster_size = (k/n_steps) * max_cluster_size;

    for j = 1 : n_steps

        subspace_dim = (j/n_steps) * max_dimension;

        truth = reshape(repmat(1:n_space, cluster_size, 1), n_space * cluster_size, 1);

        for i = 1 : n_runs

            for a = 1 : length(ambient_list)

                ambient_dim = ambient_list(a);

                base_basis = orth(randn(ambient_dim, subspace_dim / 3));

                tic;

                A = zeros(ambient_dim, n_space *cluster_size);

                for t = 1 : n_space

                    basis = [base_basis orth(randn(ambient_dim, subspace_dim - (subspace_dim/3)))];

                    A(:, cluster_size*(t-1) + 1: cluster_size*t) = basis * (rand(subspace_dim, cluster_size) - 0.5);

                end

                toc;

                X = A;

                X_normed = normalize(X);

                tic;
                [Z_ssc_relaxed] = ssc_relaxed(X_normed, 0.1);
                ssc_runtime(k, j, i, a) = toc;

                clusters_ssc = condense_clusters(ncutW((abs(Z_ssc_relaxed)+abs(Z_ssc_relaxed')), n_space),1);
                ssc_accuracy(k, j, i, a) = Misclassification(clusters_ssc, truth);

                tic;
                [Z_kssc_relaxed] = kssc_relaxed_par(X_normed, 0.1, min(round(subspace_dim * 1.5), round(cluster_size/2)));
                kssc_runtime(k, j, i, a) = toc;

                try
                clusters_kssc = condense_clusters(ncutW((abs(Z_kssc_relaxed)+abs(Z_kssc_relaxed')), n_space),1);
                kssc_accuracy(k, j, i, a) = Misclassification(clusters_kssc, truth);
                catch
                    disp('oops');
                end

            end
                
        end

    end

end

matlabpool close

save('test_synthetic_dvsn');

cbar_max = max([ssc_accuracy(:); kssc_accuracy(:)]);

for k = 1 : 3

    imagesc(flipud((mean(ssc_accuracy(:,:,:,k), 3))'))

    colorbar
    caxis([0, cbar_max]);

    set(gca, 'XTick', 1:10)
    set(gca, 'XTickLabel', 15:15:150)

    set(gca, 'YTick', 1:10)
    set(gca, 'YTickLabel', fliplr(3:3:30))

    xlabel('Cluster Size', 'FontSize', 18);

    ylabel('Subspace Dimension', 'FontSize', 18);

    set(gcf,'renderer','opengl');
    print(gcf, '-depsc2', ['test_syn_dvsn_ssc' num2str(k) '.eps']);
    
    close all
    
    imagesc(flipud((mean(kssc_accuracy(:,:,:,k), 3))'))

    colorbar
    caxis([0, cbar_max]);

    set(gca, 'XTick', 1:10)
    set(gca, 'XTickLabel', 15:15:150)

    set(gca, 'YTick', 1:10)
    set(gca, 'YTickLabel', fliplr(3:3:30))

    xlabel('Cluster Size', 'FontSize', 18);

    ylabel('Subspace Dimension', 'FontSize', 18);

    set(gcf,'renderer','opengl');
    print(gcf, '-depsc2', ['test_syn_dvsn_kssc' num2str(k) '.eps']);
    
    close all
    

end

