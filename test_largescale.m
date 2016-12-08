paths = ['common:', genpath('libs'), genpath('SSSC')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end


rng(1)

load('data/TIRlib');

n_steps = 40;
n_spectra = 3;

missrate_ssc = zeros(n_steps, 1);
missrate_kssc = zeros(n_steps, 1);
missrate_gsc = zeros(n_steps, 1);
missrate_tsc = zeros(n_steps, 1);
missrate_sssc = zeros(n_steps, 1);
missrate_gfs = zeros(n_steps, 1);

time_ssc = zeros(n_steps, 1);
time_kssc = zeros(n_steps, 1);
time_gsc = zeros(n_steps, 1);
time_tsc = zeros(n_steps, 1);
time_sssc = zeros(n_steps, 1);
time_gfs = zeros(n_steps, 1);

n_cluster = 5;

cluster_size_list = (100:100:4000)';

for cur_step = 1 : n_steps

    disp([int2str(cur_step) ' of ' int2str(n_steps)]);
    
    cluster_size = cluster_size_list(cur_step);
    
    B = zeros(321, n_cluster * cluster_size);
    
    for j = 1 : n_cluster
        spectra_indices = randi(size(A,2),n_spectra,1);

        B(:,((j-1)*cluster_size)+1:j*cluster_size) = A(:, spectra_indices) * gen_depcoeff(0.1, 0.001, n_spectra, cluster_size);
    end

    truth = reshape(repmat(1:n_cluster, cluster_size, 1), n_cluster * cluster_size, 1);
    
    X_normed = normalize(B);

    % Robust Subspace Clustering via Thresholding
    tic;
    Z_tsc = tsc_sc(X_normed, 20);
    time_tsc(cur_step, 1) = toc;

    [tsc_clusters,~,~] = ncutW((abs(Z_tsc)+abs(Z_tsc')), n_cluster);
    clusters_tsc = condense_clusters(tsc_clusters,1);
    missrate_tsc(cur_step, 1) = Misclassification(clusters_tsc, truth);

    % Scalable SSC
    tic;
    try
        clusters_sssc = get_sssc_clusters(X_normed, 20, 0.1, n_cluster );
        missrate_sssc(cur_step, 1) = Misclassification(clusters_sssc, truth);
    catch
        missrate_sssc(cur_step, 1) = 1;
    end
    time_sssc(cur_step, 1) = toc;

    % kSSC
    tic;
    [Z_kssc] = kssc_relaxed_par(X_normed, 0.1, 20);
    time_kssc(cur_step, 1) = toc;

    [kssc_clusters,~,~] = ncutW((abs(Z_kssc)+abs(Z_kssc')), n_cluster);
    clusters_kssc = condense_clusters(kssc_clusters,1);
    missrate_kssc(cur_step, 1) = Misclassification(clusters_kssc, truth);
    
    if (cluster_size <= 500)
    
        % SSC
        tic;
        [Z_ssc] = ssc_relaxed(X_normed, 0.1);
        time_ssc(cur_step, 1) = toc;

        [ssc_clusters,~,~] = ncutW((abs(Z_ssc)+abs(Z_ssc')), n_cluster);
        clusters_ssc = condense_clusters(ssc_clusters,1);
        missrate_ssc(cur_step, 1) = Misclassification(clusters_ssc, truth);
        
        % Greedy Feature Selection
        tic;
        Z_gfs = gfs_sc(X_normed, 10);
        time_gfs(cur_step, 1) = toc;

        [gfs_clusters,~,~] = ncutW((abs(Z_gfs)+abs(Z_gfs')), n_cluster);
        clusters_gfs = condense_clusters(gfs_clusters,1);
        missrate_gfs(cur_step, 1) = Misclassification(clusters_gfs, truth);
    
    end

    
    if (cluster_size <= 1000)
        
        % Greedy Subspace Clustering
        tic;
        Z_gsc = gsc_sc(X_normed, 10, 5);
        time_gsc(cur_step, 1) = toc;    

        [gsc_clusters,~,~] = ncutW((abs(Z_gsc)+abs(Z_gsc')), n_cluster);
        clusters_gsc = condense_clusters(gsc_clusters,1);
        missrate_gsc(cur_step, 1) = Misclassification(clusters_gsc, truth);
        
    end
    
    
end

save('test_largescale');

h1 = plot(time_ssc(1:5), '-*m');
hold on
h2 = plot(time_gfs(1:5), '-*b');
h3 = plot(time_gsc(1:10), '-*r');
h5 = plot(time_tsc, '-*k', 'color', [1, .5, 0]);
h4 = plot(time_sssc, '-*g');
h6 = plot(time_kssc, '-*k');

set(gca, 'fontsize', 14);

legend([h1, h2, h3, h4, h5, h6], 'SSC', 'GFS', 'GSC', 'TSC', 'SSSC', 'kSSC', 'Location','northeast')

xlim([1 40]);
% ylim([0 100]);

set(gca, 'XTick', [1, 10, 20, 30, 40])
set(gca, 'XTickLabel', 100:975:4000)

xlabel('Cluster Size', 'FontSize', 18);

ylabel('Running Time (seconds)', 'FontSize', 18);

print(gcf, '-depsc2', 'test_largescale.eps');

close all
