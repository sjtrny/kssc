paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

rng(1);

n_cores = 8;

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool(n_cores);
end

% matlabpool(n_cores);

rows = 50;
n_space = 5;
cluster_size = 20;
subspace_dim = 10;

noise_mag_list = [0, 0.05, 0.1];

n_steps = 10;
reps = 50;

ssc_accuracy = zeros(n_steps, n_steps, reps, length(noise_mag_list));
kssc_accuracy = zeros(n_steps, n_steps, reps, length(noise_mag_list));

ssc_runtime = zeros(n_steps, n_steps, reps, length(noise_mag_list));
kssc_runtime = zeros(n_steps, n_steps, reps, length(noise_mag_list));

psnr_list = zeros(n_steps, n_steps, reps, length(noise_mag_list));

max_mean = 1;
max_var = 5;

truth = reshape(repmat(1:n_space, cluster_size, 1), n_space *cluster_size, 1);

for k = 1 : n_steps
    
    m = (k/n_steps) * max_mean;
    
    for j = 1 : n_steps
        
        v = (j/n_steps) * max_var;
        
        for i = 1 : reps
            
            for a = 1 : length(noise_mag_list)
            
                tic;
                A = gen_depmultivar_data(rows, subspace_dim, cluster_size, n_space, m, v);
                toc;
                
                B = normalize(A);
                
                w = randn(size(B)) * noise_mag_list(a);
                X = B + w;
                
                psnr_list(k, j, i, a) = psnr(B, X);

                X_normed = normalize(X);

                tic;
                [Z_ssc_relaxed] = ssc_relaxed(X_normed, 0.1);
                ssc_runtime(k, j, i, a) = toc;

                clusters_ssc = condense_clusters(ncutW((abs(Z_ssc_relaxed)+abs(Z_ssc_relaxed')), n_space),1);
                ssc_accuracy(k, j, i, a) = Misclassification(clusters_ssc, truth);

                tic;
                [Z_kssc_relaxed] = kssc_relaxed_par(X_normed, 0.1, 10);
                kssc_runtime(k, j, i, a) = toc;

                clusters_kssc = condense_clusters(ncutW((abs(Z_kssc_relaxed)+abs(Z_kssc_relaxed')), n_space),1);
                kssc_accuracy(k, j, i, a) = Misclassification(clusters_kssc, truth);
            
            end

        end
        
    
    end
    
end

matlabpool close

save('test_synthetic_mv');

cbar_max = max([ssc_accuracy(:); kssc_accuracy(:)]);

for k = 1 : 3

    imagesc(flipud(mean(ssc_accuracy(:,:,:,k), 3)))

    caxis([0, cbar_max]);
    colorbar

    set(gca, 'XTick', 1:10)
    set(gca, 'XTickLabel', 0.5:0.5:5)

    set(gca, 'YTick', 1:10)
    set(gca, 'YTickLabel', fliplr(0.1:0.1:1))

    xlabel('Variance', 'FontSize', 18);

    ylabel('Mean', 'FontSize', 18);

    set(gcf,'renderer','opengl');
    print(gcf, '-depsc2', ['test_syn_mv_ssc' num2str(k) '.eps']);
    
    close all
    
    
    
    imagesc(flipud(mean(kssc_accuracy(:,:,:,k), 3)))

    caxis([0, cbar_max]);
    colorbar

    set(gca, 'XTick', 1:10)
    set(gca, 'XTickLabel', 0.5:0.5:5)

    set(gca, 'YTick', 1:10)
    set(gca, 'YTickLabel', fliplr(0.1:0.1:1))

    xlabel('Variance', 'FontSize', 18);

    ylabel('Mean', 'FontSize', 18);

    set(gcf,'renderer','opengl');
    print(gcf, '-depsc2', ['test_syn_mv_kssc' num2str(k) '.eps']);
    
    close all
    

end



