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

max_nn = 500;
n_steps = 10;

run_time_kssc = zeros(n_steps, 2);

cluster_size = 200;

tic;
A = gen_depmultivar_data(rows, 4, cluster_size, n_space, 0.25, 0.001);
toc;

w = randn(size(A)) * 0.05;
X = A + w;

X_normed = normalize(X);

nn_list = zeros(n_steps, 1);

for k = 1 : n_steps
    
    nn_list(k) = (k/n_steps) * max_nn;
    
    tic;
    [Z_kssc_exact] = kssc_exact_par(X_normed, 0.1, nn_list(k));
    run_time_kssc(k, 1) = toc;
    
    tic;
    [Z_kssc_relaxed] = kssc_relaxed_par(X_normed, 0.1, nn_list(k));
    run_time_kssc(k, 2) = toc;
    
end

matlabpool close

save('k_exact_relaxed');

% figure
% h1 = plot(run_time_ssc(:,1), '-*b');
% hold
% h2 = plot(run_time_ssc(:,2), '-ob');
% 
% legend([h1, h2], 'SSC Exact', 'SSC Relaxed', 'Location', 'NorthWest');
% 
% set(gca, 'fontsize', 14);
% 
% xlabel('Number of Processing Cores/Threads', 'FontSize', 18);
% ylabel('Running Time (seconds)', 'FontSize', 18);
% 
% print(gcf, '-depsc2', 'parcores_runtime_ssc.eps');
% 
% close all
% 
% figure
% h3 = plot(run_time_kssc(:,1), '-*r');
% hold
% h4 = plot(run_time_kssc(:,2), '-or');
% 
% legend([h3, h4], 'kSSC Exact', 'kSSC Relaxed', 'Location', 'NorthWest');
% 
% set(gca, 'fontsize', 14);
% 
% xlabel('Number of Processing Cores/Threads', 'FontSize', 18);
% ylabel('Running Time (seconds)', 'FontSize', 18);
% 
% print(gcf, '-depsc2', 'parcores_runtime_kssc.eps');
% 
% close all
