paths = ['common:', genpath('libs'), genpath('SSSC')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end

rng(1);

n_runs = 50;

n_steps = 10;
max_noise = 0.4;

missrate_ssc = zeros(1, n_runs);
missrate_kssc = zeros(1, n_runs);
missrate_gsc = zeros(1, n_runs);
missrate_tsc = zeros(1, n_runs);
missrate_sssc = zeros(1, n_runs);

time_ssc = zeros(1, n_runs);
time_kssc = zeros(1, n_runs);
time_gsc = zeros(1, n_runs);
time_tsc = zeros(1, n_runs);
time_sssc = zeros(1, n_runs);

available_class_inds = [1:13, 15:39];

n_pic = 64;
n_classes = 3;
n_points = n_pic * n_classes;

truth = reshape(repmat(1:n_classes, n_pic, 1), n_points, 1);
    
for k = 1 : n_runs

    disp([int2str(k) ' of ' int2str(n_runs)]);

    class_inds_local = randi(length(available_class_inds), n_classes, 1);
    class_inds_global = available_class_inds(class_inds_local);

%         B = zeros(32256, n_points);
    B = zeros(96*84, n_points);

    for j = 1 : n_classes

        k_list = dir(['data/CroppedYale/yaleB' num2str(class_inds_global(j), '%02.0f')  '/*P00A*.pgm*']);

        for i = 1 : n_pic

%                 B(:, j*n_pic - n_pic + i) = reshape( double( imresize(imread(['data/CroppedYale/yaleB' num2str(class_indices(j), '%02.0f') '/' k_list(i).name]), [42 48]))/255, 32256, 1);
            B(:, j*n_pic - n_pic + i) = reshape( double( imresize(imread(['data/CroppedYale/yaleB' num2str(class_inds_global(j), '%02.0f') '/' k_list(i).name]), [96 84]))/255, 96*84, 1);


        end

    end

    X_normed = normalize(B);

    % SSC
    tic;
    [Z_ssc] = ssc_relaxed(X_normed, 0.1);
    time_ssc(1, k) = toc;

    [ssc_clusters,~,~] = ncutW((abs(Z_ssc)+abs(Z_ssc')), n_classes);
    clusters_ssc = condense_clusters(ssc_clusters,1);
    missrate_ssc(1, k) = Misclassification(clusters_ssc, truth);

    % kSSC
    tic;
    [Z_kssc] = kssc_relaxed_par(X_normed, 0.1, 25);
    time_kssc(1, k) = toc;

    [kssc_clusters,~,~] = ncutW((abs(Z_kssc)+abs(Z_kssc')), n_classes);
    clusters_kssc = condense_clusters(kssc_clusters,1);
    missrate_kssc(1, k) = Misclassification(clusters_kssc, truth);

    % Greedy Subspace Clustering
    tic;
    Z_gsc = gsc_sc(X_normed, 25, 9);
    time_gsc(1, k) = toc;

    [gsc_clusters,~,~] = ncutW((abs(Z_gsc)+abs(Z_gsc')), n_classes);
    clusters_gsc = condense_clusters(gsc_clusters,1);
    missrate_gsc(1, k) = Misclassification(clusters_gsc, truth);

    % Robust Subspace Clustering via Thresholding
    tic;
    Z_tsc = tsc_sc(X_normed, 25);
    time_tsc(1, k) = toc;

    [tsc_clusters,~,~] = ncutW((abs(Z_tsc)+abs(Z_tsc')), n_classes);
    clusters_tsc = condense_clusters(tsc_clusters,1);
    missrate_tsc(1, k) = Misclassification(clusters_tsc, truth);

    % Scalable SSC
    tic;
    try
        clusters_sssc = get_sssc_clusters(X_normed, 25, 0.1, n_classes );
        missrate_sssc(1, k) = Misclassification(clusters_sssc, truth);
    catch
        missrate_sssc(1, k) = 1;
    end
    time_sssc(1, k) = toc;


end

save('test_faces');


mean_col = (mean([missrate_ssc; missrate_gsc; missrate_tsc; missrate_sssc; missrate_kssc], 2)) * 100;

median_col = (median([missrate_ssc; missrate_gsc; missrate_tsc; missrate_sssc; missrate_kssc], 2)) * 100;

min_col = (min([missrate_ssc; missrate_gsc; missrate_tsc; missrate_sssc; missrate_kssc], [], 2)) * 100;

max_col = (max([missrate_ssc; missrate_gsc; missrate_tsc; missrate_sssc; missrate_kssc], [], 2)) * 100;

std_col = std([missrate_ssc; missrate_gsc; missrate_tsc; missrate_sssc; missrate_kssc], [], 2) * 100;

mean_time = mean([time_ssc; time_kssc; time_tsc; time_sssc; time_kssc], 2);

names = {'SSC', 'GSC', 'TSC', 'SSSC', 'kSSC'};


for table_cnt = 1 : 5

    fprintf('%s & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.1f\\%% & %.2f \\\\ \n', names{table_cnt}, mean_col(table_cnt), median_col(table_cnt), min_col(table_cnt), max_col(table_cnt), std_col(table_cnt), mean_time(table_cnt));
    
end

