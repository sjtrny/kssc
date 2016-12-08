paths = ['common:', genpath('libs')];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

rng(1);

max_cores = 8;

rows = 500;
n_space = 5;
cluster_size = 200;

tic;
A = gen_depmultivar_data(rows, 4, cluster_size, n_space, 0.25, 0.001);
toc;

w = randn(size(A)) * 0.05;
X = A + w;

X_normed = normalize(X);

run_time_ssc = zeros(max_cores, 2);
run_time_kssc = zeros(max_cores, 2);

for k = 1 : max_cores
    
    matlabpool(k);
    
    tic;
    [Z_ssc_exact] = ssc_exact_par(X_normed, 0.1);
    run_time_ssc(k, 1) = toc;
    
    tic;
    [Z_ssc_relaxed] = ssc_relaxed_par(X_normed, 0.1);
    run_time_ssc(k, 2) = toc;
    
    tic;
    [Z_kssc_exact] = kssc_exact_par(X_normed, 0.1, 6);
    run_time_kssc(k, 1) = toc;
    
    tic;
    [Z_kssc_relaxed] = kssc_relaxed_par(X_normed, 0.1, 6);
    run_time_kssc(k, 2) = toc;
   
    matlabpool close
    
end

save('parcores_runtime');

figure
h1 = plot(run_time_ssc(:,1), '-ob');
hold
h2 = plot(run_time_ssc(:,2), '-*g');
h3 = plot(run_time_kssc(:,1), '-*r');
h4 = plot(run_time_kssc(:,2), '-ok');

legend([h1, h2, h3, h4], 'SSC Exact (Column Parallel)', 'SSC Relaxed (Column Parallel)', 'kSSC Exact (Column Parallel)', 'kSSC Relaxed (Column Parallel)', 'Location', 'NorthWest');

set(gca, 'fontsize', 14);

ylim([-10, 300]);

xlabel('Number of Processing Cores/Threads', 'FontSize', 18);
ylabel('Running Time (seconds)', 'FontSize', 18);

print(gcf, '-depsc2', 'parcoresruntimeoverall.eps');

close all

figure
h5 = plot(run_time_kssc(:,1), '-*r');
hold
h6 = plot(run_time_kssc(:,2), '-ok');

legend([h5, h6], 'kSSC Exact (Column Parallel)', 'kSSC Relaxed (Column Parallel)', 'Location', 'NorthEast');

set(gca, 'fontsize', 14);

ylim([0, 1.5]);

xlabel('Number of Processing Cores/Threads', 'FontSize', 18);
ylabel('Running Time (seconds)', 'FontSize', 18);

print(gcf, '-depsc2', 'parcoresruntimekssc.eps');

close all
