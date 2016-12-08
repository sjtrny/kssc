paths = ['common:', genpath('libs'), genpath('SSSC')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end

rng(1)

file_list = dir('data/Hopkins155/*.mat');

n_captures = length(file_list);

n_steps = 10;
max_noise = 0.4;

missrate_ssc = zeros(n_steps + 1, n_captures);
missrate_kssc = zeros(n_steps + 1, n_captures);
missrate_gfs = zeros(n_steps + 1, n_captures);
missrate_gsc = zeros(n_steps + 1, n_captures);
missrate_tsc = zeros(n_steps + 1, n_captures);
missrate_sssc = zeros(n_steps + 1, n_captures);

time_ssc = zeros(n_steps + 1, n_captures);
time_kssc = zeros(n_steps + 1, n_captures);
time_gfs = zeros(n_steps + 1, n_captures);
time_gsc = zeros(n_steps + 1, n_captures);
time_tsc = zeros(n_steps + 1, n_captures);
time_sssc = zeros(n_steps + 1, n_captures);


psnr_list = zeros(n_steps + 1, n_captures);

for noise_step = 1 : n_steps + 1

    noise_mag = max(((noise_step-1)/n_steps)*max_noise, 0);
    
    for k = 1 : n_captures

        disp([int2str(k) ' of ' int2str(n_captures)]);
        
        load(['data/Hopkins155/' file_list(k).name]);

        n_cluster = length(unique(s));
    %     n_samples = points;

        X = zeros(2*frames, points);
        for f = 1 : frames
            X(f*2 - 1 : f*2, :) = x(1:2, :, f);
        end

        [a, b] = sort(s);

        truth = s(b);

        X = X(:,b);

        noise = randn(size(X));

        w = noise * noise_mag;

        psnr_list(noise_step, k) = psnr(X, X + w);
        
        X_normed = normalize(X + w);

        % SSC
        tic;
        [Z_ssc] = ssc_relaxed(X_normed, 0.1);
        time_ssc(noise_step, k) = toc;
    
        [ssc_clusters,~,~] = ncutW((abs(Z_ssc)+abs(Z_ssc')), n_cluster);
        clusters_ssc = condense_clusters(ssc_clusters,1);
        missrate_ssc(noise_step, k) = Misclassification(clusters_ssc, truth);
    
        % Greedy Feature Selection
        tic;
        Z_gfs = gfs_sc(X_normed, 10);
        time_gfs(noise_step, k) = toc;
    
        [gfs_clusters,~,~] = ncutW((abs(Z_gfs)+abs(Z_gfs')), n_cluster);
        clusters_gfs = condense_clusters(gfs_clusters,1);
        missrate_gfs(noise_step, k) = Misclassification(clusters_gfs, truth);
    
        % Greedy Subspace Clustering
        tic;
        Z_gsc = gsc_sc(X_normed, 10, 5);
        time_gsc(noise_step, k) = toc;
    
        [gsc_clusters,~,~] = ncutW((abs(Z_gsc)+abs(Z_gsc')), n_cluster);
        clusters_gsc = condense_clusters(gsc_clusters,1);
        missrate_gsc(noise_step, k) = Misclassification(clusters_gsc, truth);

        % Robust Subspace Clustering via Thresholding
        tic;
        Z_tsc = tsc_sc(X_normed, 20);
        time_tsc(noise_step, k) = toc;

        [tsc_clusters,~,~] = ncutW((abs(Z_tsc)+abs(Z_tsc')), n_cluster);
        clusters_tsc = condense_clusters(tsc_clusters,1);
        missrate_tsc(noise_step, k) = Misclassification(clusters_tsc, truth);
        
        % Scalable SSC
        tic;
        try
            clusters_sssc = get_sssc_clusters(X_normed, 20, 0.1, n_cluster );
            missrate_sssc(noise_step, k) = Misclassification(clusters_sssc, truth);
        catch
            missrate_sssc(noise_step, k) = 1;
        end
        time_sssc(noise_step, k) = toc;
        
        % kSSC
        tic;
        [Z_kssc] = kssc_relaxed_par(X_normed, 0.1, 20);
        time_kssc(noise_step, k) = toc;

        [kssc_clusters,~,~] = ncutW((abs(Z_kssc)+abs(Z_kssc')), n_cluster);
        clusters_kssc = condense_clusters(kssc_clusters,1);
        missrate_kssc(noise_step, k) = Misclassification(clusters_kssc, truth);


    end



end

save('test_motion');


mean_error = [mean(missrate_ssc, 2), mean(missrate_gfs, 2), mean(missrate_gsc, 2), mean(missrate_tsc, 2), mean(missrate_sssc, 2), mean(missrate_kssc, 2)];
median_error = [median(missrate_ssc, 2), median(missrate_gfs, 2), median(missrate_gsc, 2), median(missrate_tsc, 2), median(missrate_sssc, 2), median(missrate_kssc, 2)];
max_error =  [max(missrate_ssc, [], 2), max(missrate_gfs, [], 2), max(missrate_gsc, [], 2), max(missrate_tsc, [], 2), max(missrate_sssc, [], 2), max(missrate_kssc, [], 2)];
min_error =  [min(missrate_ssc, [], 2), min(missrate_gfs, [], 2), min(missrate_gsc, [], 2), min(missrate_tsc, [], 2), min(missrate_sssc, [], 2), min(missrate_kssc, [], 2)];

median_time = [median(time_ssc, 2), median(time_gfs, 2), median(time_gsc, 2), median(time_tsc, 2), median(time_sssc, 2), median(time_kssc, 2)];

mean_psnr = mean(psnr_list, 2);

mean_error = mean_error * 100;
median_error = median_error * 100;
max_error = max_error * 100;
min_error = min_error * 100;

% Mean
h1 = plot(1:11, mean_error(:,1), '-*m');
hold on
h2 = plot(1:11, mean_error(:,2), '-*b');
h3 = plot(1:11, mean_error(:,3), '-*r');
h4 = plot(1:11, mean_error(:,4), '-*k', 'color', [1, .5, 0]);
h5 = plot(1:11, mean_error(:,5), '-*g');
h6 = plot(1:11, mean_error(:,6), '-*k');

set(gca, 'fontsize', 14);

legend([h1, h2, h3, h4, h5, h6], 'SSC', 'GFS', 'GSC', 'TSC', 'SSSC', 'kSSC', 'Location','northwest')

xlim([1 11]);
ylim([0 100]);

set(gca, 'XTick', 1:11)
set(gca, 'XTickLabel', mean(round(mean_psnr), 2))

xlabel('PSNR', 'FontSize', 18);

ylabel('Subspace Clustering Error (SCE)', 'FontSize', 18);

print(gcf, '-depsc2', 'motion_mean.eps');

close all

% Median
h1 = plot(1:11, median_error(:,1), '-*m');
hold on
h2 = plot(1:11, median_error(:,2), '-*b');
h3 = plot(1:11, median_error(:,3), '-*r');
h4 = plot(1:11, median_error(:,4), '-*k', 'color', [1, .5, 0]);
h5 = plot(1:11, median_error(:,5), '-*g');
h6 = plot(1:11, median_error(:,6), '-*k');

set(gca, 'fontsize', 14);

legend([h1, h2, h3, h4, h5, h6], 'SSC', 'GFS', 'GSC', 'TSC', 'SSSC', 'kSSC', 'Location','northwest')

xlim([1 11]);
ylim([0 100]);

set(gca, 'XTick', 1:11)
set(gca, 'XTickLabel', mean(round(mean_psnr), 2))

xlabel('PSNR', 'FontSize', 18);

ylabel('Subspace Clustering Error (SCE)', 'FontSize', 18);

print(gcf, '-depsc2', 'motion_median.eps');

close all

% Max
h1 = plot(1:11, max_error(:,1), '-*m');
hold on
h2 = plot(1:11, max_error(:,2), '-*b');
h3 = plot(1:11, max_error(:,3), '-*r');
h4 = plot(1:11, max_error(:,4), '-*k', 'color', [1, .5, 0]);
h5 = plot(1:11, max_error(:,5), '-*g');
h6 = plot(1:11, max_error(:,6), '-*k');

set(gca, 'fontsize', 14);

legend([h1, h2, h3, h4, h5, h6], 'SSC', 'GFS', 'GSC', 'TSC', 'SSSC', 'kSSC', 'Location','northwest')

xlim([1 11]);
ylim([0 100]);

set(gca, 'XTick', 1:11)
set(gca, 'XTickLabel', mean(round(mean_psnr), 2))

xlabel('PSNR', 'FontSize', 18);

ylabel('Subspace Clustering Error (SCE)', 'FontSize', 18);

print(gcf, '-depsc2', 'motion_max.eps');

close all

% Min
h1 = plot(1:11, min_error(:,1), '-*m');
hold on
h2 = plot(1:11, min_error(:,2), '-*b');
h3 = plot(1:11, min_error(:,3), '-*r');
h4 = plot(1:11, min_error(:,4), '-*k', 'color', [1, .5, 0]);
h5 = plot(1:11, min_error(:,5), '-*g');
h6 = plot(1:11, min_error(:,6), '-*k');

set(gca, 'fontsize', 14);

legend([h1, h2, h3, h4, h5, h6], 'SSC', 'GFS', 'GSC', 'TSC', 'SSSC', 'kSSC', 'Location','northwest')

xlim([1 11]);
ylim([0 100]);

set(gca, 'XTick', 1:11)
set(gca, 'XTickLabel', mean(round(mean_psnr), 2))

xlabel('PSNR', 'FontSize', 18);

ylabel('Subspace Clustering Error (SCE)', 'FontSize', 18);

print(gcf, '-depsc2', 'motion_min.eps');

close all

% Running time

h1 = plot(1:11, median_time(:,1), '-*m');
hold on
h2 = plot(1:11, median_time(:,2), '-*b');
h3 = plot(1:11, median_time(:,3), '-*r');
h4 = plot(1:11, median_time(:,4), '-*k', 'color', [1, .5, 0]);
h5 = plot(1:11, median_time(:,5), '-*g');
h6 = plot(1:11, median_time(:,6), '-*k');

set(gca, 'fontsize', 14);

legend([h1, h2, h3, h4, h5, h6], 'SSC', 'GFS', 'GSC', 'TSC', 'SSSC', 'kSSC', 'Location','northwest')

xlim([1 11]);
ylim([0 4]);

set(gca, 'XTick', 1:11)
set(gca, 'XTickLabel', mean(round(mean_psnr), 2))

xlabel('PSNR', 'FontSize', 18);

ylabel('Running Time (seconds)', 'FontSize', 18);

print(gcf, '-depsc2', 'motion_time.eps');

close all

