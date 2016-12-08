paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

rng(1);

rows = 500;
n_space = 5;

n_cores = 8;

if ~exist('matlabpool','file') && ~matlabpool('size') == 0
    matlabpool close
end

matlabpool(n_cores);

max_cluster_size = 200;
n_steps = 10;

run_time_ssc = zeros(n_steps, 4);
run_time_kssc = zeros(n_steps, 2);

for k = 1 : n_steps
    
    cluster_size = (k/n_steps)* max_cluster_size;
    
    tic;
    A = gen_depmultivar_data(rows, 4, cluster_size, n_space, 0.25, 0.001);
    toc;
    
    w = randn(size(A)) * 0.05;
    X = A + w;
    
    X_normed = normalize(X);
    
    tic;
    [Z_ssc_exact] = ssc_exact(X_normed, 0.1);
    run_time_ssc(k, 1) = toc;
    
    tic;
    [Z_ssc_relaxed] = ssc_relaxed(X_normed, 0.1);
    run_time_ssc(k, 2) = toc;
    
    tic;
    [Z_ssc_exact_par] = ssc_exact_par(X_normed, 0.1);
    run_time_ssc(k, 3) = toc;
    
    tic;
    [Z_ssc_relaxed_par] = ssc_relaxed_par(X_normed, 0.1);
    run_time_ssc(k, 4) = toc;
    
    tic;
    [Z_kssc_exact] = kssc_exact_par(X_normed, 0.1, 6);
    run_time_kssc(k, 1) = toc;
    
    tic;
    [Z_kssc_relaxed] = kssc_relaxed_par(X_normed, 0.1, 6);
    run_time_kssc(k, 2) = toc;
    
end

matlabpool close

save('runtime');

figure
h1 = plot(run_time_ssc(:,1), '-*b');
hold
h2 = plot(run_time_ssc(:,3), '-ob');
h3 = plot(run_time_ssc(:,2), '-*g');
h4 = plot(run_time_ssc(:,4), '-og');

h5 = plot(run_time_kssc(:,1), '-*r');
h6 = plot(run_time_kssc(:,2), '-ok');

ylim([-10, 250]);
xlim([1, 10]);

legend([h1, h2, h3, h4, h5, h6], 'SSC Exact', 'SSC Exact (Column Parallel)', 'SSC Relaxed', 'SSC Relaxed (Column Parallel)', 'kSSC Exact (Column Parallel)', 'kSSC Relaxed (Column Parallel)', 'Location', 'NorthWest');

set(gca, 'fontsize', 14);

set(gca, 'XTickLabel', (20:20:200)*5);

xlabel('N - Number of datapoints', 'FontSize', 18);
ylabel('Running Time (seconds)', 'FontSize', 18);

print(gcf, '-depsc2', 'runtimeoverall.eps');

close all

figure
h7 = plot(run_time_kssc(:,1), '-*r');
hold
h8 = plot(run_time_kssc(:,2), '-ok');

ylim([0, 0.8]);
xlim([1, 10]);

legend([h7, h8], 'kSSC Exact (Column Parallel)', 'kSSC Relaxed (Column Parallel)', 'Location', 'NorthWest');

set(gca, 'fontsize', 14);

set(gca, 'XTickLabel', (20:20:200)*5);

xlabel('N - Number of datapoints', 'FontSize', 18);
ylabel('Running Time (seconds)', 'FontSize', 18);

print(gcf, '-depsc2', 'runtimekssc.eps');

close all
